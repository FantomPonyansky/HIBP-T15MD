'''
Heavy Ion Beam Probe partile tracing library
'''
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
import hibpcalc.geomfunc as gf
from hibpcalc.fields import return_E
#from hibpcalc.fields import _return_E as return_E

from hibpcalc.misc import runge_kutt, argfind_rv, find_fork


# %% define class for trajectories
class Traj():
    '''
    Trajectory object
    '''

    def __init__(self, q, m, Ebeam, r0, alpha, beta, U, dt=1e-7):
        '''

        Parameters
        ----------
        q : float
            particle charge [Co]
        m : float
            particle mass [kg]
        Ebeam : float
            beam energy [keV]
        r0 : np.array
            initial point of the trajectory [m]
        alpha : float
            injection angle in XY plane [rad]
        beta : float
            injection angle in XZ plane [rad]
        U : dict
            dict of voltages in [kV] keys=[A1 B1 A2 B2 A3 B3 an]
        dt : float, optional
            timestep for RK algorithm [s]. The default is 1e-7.

        Returns
        -------
        None.

        '''
        self.q = q
        self.m = m
        self.Ebeam = Ebeam
        # particle velocity:
        Vabs = np.sqrt(2 * Ebeam * 1.602176634E-16 / m)
        V0 = gf.calc_vector(-Vabs, alpha, beta)
        self.alpha = alpha
        self.beta = beta
        self.U = U
        self.RV0 = np.array([np.hstack((r0, V0))])  # initial condition
        # array with r,V for the primary trajectory
        self.RV_prim = self.RV0
        self.tag_prim = [1]
        # array with r,V for the secondary trajectory
        self.RV_sec = np.array([[]])
        self.tag_sec = [2]
        # list to contain RV of the whole fan:
        self.Fan = []
        # time step for primary orbit:
        self.dt1 = dt
        # time step for secondary orbit:
        self.dt2 = dt
        # flags
        self.IsAimXY = False
        self.IsAimZ = False
        self.fan_ok = False
        self.IntersectGeometry = {'A2': False, 'B2': False, 'chamb': False}
        self.IntersectGeometrySec = {'A3': False, 'B3': False, 'A4': False,
                                     'chamb': False}
        self.B_out_of_bounds = False
        self._B_out_of_bounds_debug_info = None
        
        # multislit:
        self.fan_to_slits = []
        self.RV_sec_toslits = []
        self.ion_zones = []
        
        self.log = []

    

    def print_log(self, s): 
        self.log.append(s)
        print(s)

    def pass_prim(self, E_interp, B_interp, geom, tmax=1e-5, 
                  invisible_wall_x=5.5):
        '''
        passing primary trajectory from initial point self.RV0
        E_interp : dictionary with E field interpolants
        B_interp : list with B fied interpolants
        geom : Geometry object
        '''
        print('\n Passing primary trajectory')
        # reset intersection flags
        for key in self.IntersectGeometry.keys():
            self.IntersectGeometry[key] = False
        t = 0.
        dt = self.dt1
        RV_old = self.RV0  # initial position
        RV = self.RV0  # array to collect all r, V
        k = self.q / self.m
        tag_column = [10]

        while t <= tmax:
            r = RV_old[0, :3]
            
            # Electric field:
            E_local = return_E(r, E_interp, self.U, geom)
            
            # Magnetic field:
            B_local = B_interp(r)
            if self.check_B_is_NaN(B_local, (r, 'prim')):
                break
                
            # runge-kutta step:
            RV_new = runge_kutt(k, RV_old, dt, E_local, B_local)

            # check if out of bounds for passing to aim
            if RV_new[0, 0] > invisible_wall_x and RV_new[0, 1] < 1.2:
                self.print_log('primary hit invisible wall, r = %s' % str(np.round(r, 3)))
                break
            
            # check if intersected chamber entrance
            if geom.check_chamb_intersect('prim', RV_old[0, 0:3],
                                          RV_new[0, 0:3]):
                self.print_log('Primary intersected chamber entrance')
                self.IntersectGeometry['chamb'] = True
                break

            plts_flag, plts_name = geom.check_plates_intersect(RV_old[0, 0:3],
                                                               RV_new[0, 0:3])
            if plts_flag:
                #self.print_log('Primary intersected ' + plts_name + ' plates')
                self.IntersectGeometry[plts_name] = True
                break

            if geom.check_fw_intersect(RV_old[0, 0:3], RV_new[0, 0:3]):
                self.print_log('Primary intersected first wall')
                break  # stop primary trajectory calculation
                
            # save results
            RV = np.vstack((RV, RV_new))
            tag_column = np.hstack((tag_column, 10))
            RV_old = RV_new
            t = t + dt

        else: 
            self.print_log('t <= tmax, t=%f' % t)
            # <krokhalev> debug plot
            plt.figure(347)
            self.plot()

        self.RV_prim = RV
        self.tag_prim = tag_column

    def pass_sec(self, RV0, r_aim, E_interp, B_interp, geom,
                 stop_plane_n=np.array([1., 0., 0.]), tmax=5e-5,
                 eps_xy=1e-3, eps_z=1e-3, invisible_wall_x=5.5):
        '''
        passing secondary trajectory from initial point RV0 to point r_aim
        with accuracy eps
        RV0 : initial position and velocity
        '''
        # print('Passing secondary trajectory')
        self.IsAimXY = False
        self.IsAimZ = False
        self.B_out_of_bounds = False
        self._B_out_of_bounds_debug_info = None
        # reset intersection flags
        for key in self.IntersectGeometrySec.keys():
            self.IntersectGeometrySec[key] = False
        t = 0.
        dt = self.dt2
        RV_old = RV0  # initial position
        RV = RV0  # array to collect all [r,V]
        k = 2*self.q / self.m
        tag_column = [20]

        while t <= tmax:  # to witness curls initiate tmax as 1 sec
            r = RV_old[0, :3]
            
            # Electric field:
            E_local = return_E(r, E_interp, self.U, geom)
            
            # Magnetic field:
            B_local = B_interp(r)
            if self.check_B_is_NaN(B_local, (r, 'sec')):
                break
            
            # runge-kutta step:
            RV_new = runge_kutt(k, RV_old, dt, E_local, B_local)

        
            # #check if out of bounds for passing to aim
            if RV_new[0, 0] > invisible_wall_x:
                self.print_log(
                    'secondary hit invisible wall, r = %s, invisible_wall_x = %s, r_aim = %s, stop_n = %s' % (    
                        str(np.round(r, 3)), 
                        str(np.round(invisible_wall_x, 3)), 
                        str(np.round(r_aim, 3)), 
                        str(np.round(stop_plane_n, 3))
                    ))
                plt.figure(553)
                plt.plot(RV[:, 0], RV[:, 1])
                plt.plot(RV[:, 0], RV[:, 2])

                plt.plot(self.RV_prim[:, 0], self.RV_prim[:, 1])
                plt.plot(self.RV_prim[:, 0], self.RV_prim[:, 2])
                
                # raise Exception("???")
                
                break
                
            if geom.check_chamb_intersect('sec', RV_old[0, 0:3],
                                          RV_new[0, 0:3]):
                # print('Secondary intersected chamber exit')
                self.IntersectGeometrySec['chamb'] = True

            plts_flag, plts_name = geom.check_plates_intersect(RV_old[0, 0:3], RV_new[0, 0:3])
            if plts_flag:
                self.print_log('Secondary intersected ' + plts_name + ' plates')
                self.IntersectGeometrySec[plts_name] = True

            # find last point of the secondary trajectory
            #if (RV_new[0, 0] > 2.45) and (RV_new[0, 1] < 1.5):
            if (RV_new[0, 0] > r_aim[0] - 0.085) and (RV_new[0, 1] < 1.5):
                # intersection with the stop plane:
                
                #!!! original code    
                r_intersect = gf.line_plane_intersect(stop_plane_n, r_aim, RV_new[0, :3]-RV_old[0, :3], RV_new[0, :3])
                # _r_intersect = gf.plane_segment_intersect(stop_plane_n, r_aim, RV_old[0, :3], RV_new[0, :3])
                
                
                # r_intersect = plane_segment_intersect(stop_plane_n, r_aim,
                #                                       RV_old[0, :3], RV_new[0, :3])
                # check if r_intersect is between RV_old and RV_new:
                
                #!!! original code
                if gf.is_between(RV_old[0, :3], RV_new[0, :3], r_intersect):
                
                # if not np.isnan(r_intersect):  # if r_sect is not None: 
                # if r_intersect is not None: 
                    RV_new[0, :3] = r_intersect
                    RV = np.vstack((RV, RV_new))
                    # check XY plane:
                    if (np.linalg.norm(RV_new[0, :2] - r_aim[:2]) <= eps_xy):
                        # print('aim XY!')
                        self.IsAimXY = True
                    # check XZ plane:
                    if (np.linalg.norm(RV_new[0, [0, 2]] - r_aim[[0, 2]]) <= eps_z):
                        # print('aim Z!')
                        self.IsAimZ = True
                    break
                # else: 
                #     if _r_intersect is not None: 
                #         print('WTF??')
                #         raise Exception("WTF???")
                
            # continue trajectory calculation:
            RV_old = RV_new
            t = t + dt
            RV = np.vstack((RV, RV_new))
            tag_column = np.hstack((tag_column, 20))
            # print('t secondary = ', t)

        else:
            self.print_log("max time exceeded during passing secondary")
            # <krokhalev> debug plot
            plt.figure(348)
            self.plot()
            
        self.RV_sec = RV
        self.tag_sec = tag_column

    def pass_fan(self, r_aim, E_interp, B_interp, geom,
                 stop_plane_n=np.array([1., 0., 0.]), eps_xy=1e-3, eps_z=1e-3,
                 no_intersect=False, no_out_of_bounds=False, 
                 invisible_wall_x=5.5):
        '''
        passing fan from initial point self.RV0
        '''
        print('\n Passing fan of trajectories')
        self.pass_prim(E_interp, B_interp, geom, 
                       invisible_wall_x=invisible_wall_x)
        # create a list fro secondary trajectories:
        list_sec = []
        # check intersection of primary trajectory:
        if True in self.IntersectGeometry.values():
            print('Fan list is empty')
            self.Fan = list_sec
            return

        # check eliptical radius of particle:
        # 1.5 m - major radius of a torus, elon - size along Y
        mask = np.sqrt((self.RV_prim[:, 0] - geom.R)**2 +
                       (self.RV_prim[:, 1] / geom.elon)**2) <= geom.r_plasma
        
        self.tag_prim[mask] = 11

        # list of initial points of secondary trajectories:
        # RV0_sec = self.RV_prim[(self.tag_prim == 11)]
        RV0_sec = self.RV_prim[mask]

        for RV02 in RV0_sec:
            RV02 = np.array([RV02])
            self.pass_sec(RV02, r_aim, E_interp, B_interp, geom,
                          stop_plane_n=stop_plane_n,
                          eps_xy=eps_xy, eps_z=eps_z, 
                          invisible_wall_x=invisible_wall_x)
            if (no_intersect and True in self.IntersectGeometrySec.values()) or \
               (no_out_of_bounds and self.B_out_of_bounds):
                continue
            list_sec.append(self.RV_sec)

        self.Fan = list_sec

    def pass_to_target(self, r_aim, E_interp, B_interp, geom,
                       stop_plane_n=np.array([1., 0., 0.]),
                       eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10,
                       no_intersect=False, no_out_of_bounds=False, 
                       invisible_wall_x=5.5):
        '''
        find secondary trajectory which goes directly to target
        '''
        
        if True in self.IntersectGeometry.values():
            print('There is intersection at primary trajectory')
            return
        if len(self.Fan) == 0:
            print('NO secondary trajectories')
            return
        
        # reset flags in order to let the algorithm work properly
        self.IsAimXY = False
        self.IsAimZ = False
        
        # reset intersection flags for secondaries
        for key in self.IntersectGeometrySec.keys():
            self.IntersectGeometrySec[key] = False

        # find which secondaries are higher/lower than r_aim
        # sign = -1 means higher, 1 means lower
        signs = np.array([np.sign(np.cross(RV[-1, :3], r_aim)[-1])
                          for RV in self.Fan])
        are_higher = np.argwhere(signs == -1)
        are_lower = np.argwhere(signs == 1)
        twisted_fan = False  # flag to detect twist of the fan

        if are_higher.shape[0] == 0:
            print('all secondaries are lower than aim!')
            n = int(are_lower[are_lower.shape[0]//2])
        elif are_lower.shape[0] == 0:
            print('all secondaries are higher than aim!')
            n = int(are_higher[are_higher.shape[0]//2])
        else:
            if are_higher[-1] > are_lower[0]:
                print('Fan is twisted!')
                twisted_fan = True
                n = int(are_lower[-1])
            else:
                n = int(are_higher[-1])  # find the last one which is higher
                self.fan_ok = True
        RV_old = np.array([self.Fan[n][0]])

        # find secondary, which goes directly into r_aim
        self.dt1 = self.dt1/2.
        while True:
            # make a small step along primary trajectory
            r = RV_old[0, :3]
            
            # fields
            E_local = np.array([0., 0., 0.])
            B_local = B_interp(r)
            if np.isnan(B_local).any():
                # self.print_log('Btor is nan, r = %s' % str(r))
                break
            
            #runge-kutta step
            RV_new = runge_kutt(self.q / self.m, RV_old, self.dt1,
                                E_local, B_local)
            
            # check if RV_new is in plasma
            if not (np.sqrt((RV_new[:, 0] - geom.R)**2 +
                           (RV_new[:, 1] / geom.elon)**2) <= geom.r_plasma):
                self.print_log("out of plasma during passing to target")
                break
            
            # pass new secondary trajectory
            self.pass_sec(RV_new, r_aim, E_interp, B_interp, geom,
                          stop_plane_n=stop_plane_n,
                          eps_xy=eps_xy, eps_z=eps_z, 
                          invisible_wall_x=invisible_wall_x)
            # check XY flag
            if self.IsAimXY:
                # insert RV_new into primary traj
                # find the index of the point in primary traj closest to RV_new
                ind = np.nanargmin(np.linalg.norm(self.RV_prim[:, :3] -
                                                  RV_new[0, :3], axis=1))
                if gf.is_between(self.RV_prim[ind, :3],
                              self.RV_prim[ind+1, :3], RV_new[0, :3]):
                    i2insert = ind+1
                else:
                    i2insert = ind
                self.RV_prim = np.insert(self.RV_prim, i2insert, RV_new, axis=0)
                self.tag_prim = np.insert(self.tag_prim, i2insert, 11, axis=0)
                break
            
            # check if the new secondary traj is lower than r_aim
            if (not twisted_fan and
                    np.sign(np.cross(self.RV_sec[-1, :3], r_aim)[-1]) > 0):
                # if lower, halve the timestep and try once more
                self.dt1 = self.dt1/2.
                print('dt1={}'.format(self.dt1))
                if self.dt1 < dt_min:
                    print('dt too small')
                    break
            else:
                # if higher, continue steps along the primary
                RV_old = RV_new

    def add_slits(self, n_slits):
        '''
        create empty list for secondary trajectories,
        which go to different slits
        '''
        if len(self.RV_sec_toslits) == n_slits:
            pass
        else:
            self.RV_sec_toslits = [None]*n_slits
            self.ion_zones = [None]*n_slits

    def plot_prim(self, ax, axes='XY', color='k', full_primary=False):
        '''
        plot primary trajectory
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        index = -1
        
        if min(self.RV_sec.shape) == 0:
            full_primary = True        
    
        if not full_primary:
            # find where secondary trajectory starts:
            for i in range(self.RV_prim.shape[0]):
                if np.linalg.norm(self.RV_prim[i, :3] - self.RV_sec[0, :3]) < 1e-4:
                    index = i+1
        ax.plot(self.RV_prim[:index, index_X],
                self.RV_prim[:index, index_Y],
                color=color, linewidth=2)

    def plot_sec(self, ax, axes='XY', color='r'):
        '''
        plot secondary trajectory
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        ax.plot(self.RV_sec[:, index_X], self.RV_sec[:, index_Y],
                color=color, linewidth=2)

    def plot_fan(self, ax, axes='XY', color='r'):
        '''
        plot fan of secondary trajectories
        '''
        axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
        index_X, index_Y = axes_dict[axes]
        for i in self.Fan:
            ax.plot(i[:, index_X], i[:, index_Y], color=color)
            
#%% new traj functions !!! TO BE TESTED !!!
    # def _plot_fan(self, ax, axes='XY', color='r', indexces=):
    #     '''
    #     plot fan of secondary trajectories
    #     '''
    #     axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
    #     index_X, index_Y = axes_dict[axes]
    #     for i in self.Fan:
    #         ax.plot(i[:, index_X], i[:, index_Y], color=color)   
    
    def plot(self, ax=None, axes='XY', color='r'):
        if ax is None:
            ax = plt.gca()
        self.plot_prim(ax, axes=axes, color='black')
        try:
            self.plot_sec(ax, axes=axes, color=color)
        except:
            pass
        
    def reset_sec_flags(self):
        '''
        Sets aim, b_out_of_bounds and intersection flags to False

        Returns
        -------
        None.

        '''
        self.IsAimXY = False
        self.IsAimZ = False
        self.B_out_of_bounds = False
        self._B_out_of_bounds_debug_info = None
        
        for key in self.IntersectGeometrySec.keys():
            self.IntersectGeometrySec[key] = False
            
    def check_sec_intersection(self, RV_old, RV_new):
        '''
        Checks if last step of econdadry trajectory intersect geometry or plates
        
        
        Parameters
        ----------
        RV_old : TYPE
            DESCRIPTION.
        RV_new : TYPE
            DESCRIPTION.

        Returns
        -------
        bool
            True if RV_new out of bounds
        '''

    def _pass_sec(self, RV0, r_aim, E_interp, B_interp, geom,
                 stop_plane_n=np.array([1., 0., 0.]), tmax=5e-5,
                 eps_xy=1e-3, eps_z=1e-3, invisible_wall_x=5.5, break_at_intersection=False):
        '''
        passing secondary trajectory from initial point RV0 to point r_aim
        with accuracy eps.
        #!!!
        Parameters
        ----------
        RV0 : TYPE
            Initial position and velocity
        r_aim : TYPE
            DESCRIPTION.
        E_interp : TYPE
            DESCRIPTION.
        B_interp : TYPE
            DESCRIPTION.
        geom : TYPE
            DESCRIPTION.
        stop_plane_n : np.array with shape (3,), optional
            Normal vector to the plane there secolndary trajectory will be cut.
            The default is np.array([1., 0., 0.]).
        tmax : float64, optional
            DESCRIPTION. The default is 5e-5.
        eps_xy : float64, optional
            DESCRIPTION. The default is 1e-3.
        eps_z : float64, optional
            DESCRIPTION. The default is 1e-3.
        invisible_wall_x : float64, optional
            DESCRIPTION. The default is 5.5.

        Returns
        -------
        None.

        '''
        self.reset_sec_flags()
        
        #set initial parameters        
        t = 0.
        dt = self.dt2
        RV_old = RV0  # initial position
        RV = RV0  # array to collect all [r,V]
        k = 2*self.q / self.m
        tag_column = [20]
        
        # pass trajectory
        while t <= tmax:  # to witness curls initiate tmax as 1 sec
            r = RV_old[0, :3]
            
            # fields
            E_local = return_E(r, E_interp, self.U, geom)
            B_local = B_interp(r)
            if self.check_B_is_NaN(B_local, (r, '_sec')):
                break
                
            # runge-kutta step:
            RV_new = runge_kutt(k, RV_old, dt, E_local, B_local)

            # check if out of bounds for passing to aim
            if RV_new[0, 0] > invisible_wall_x:
                self.print_log('secondary hit invisible wall, r = %s' % str(np.round(r, 3)))
                break
            
            #check if intersected geometry or plates
            if geom.check_chamb_intersect('sec', RV_old[0, 0:3], RV_new[0, 0:3]):
                self.IntersectGeometrySec['chamb'] = True
                if break_at_intersection: 
                    break

            plts_flag, plts_name = geom.check_plates_intersect(RV_old[0, 0:3],
                                                               RV_new[0, 0:3])
            if plts_flag:
                self.print_log('Secondary intersected ' + plts_name + ' plates')
                self.IntersectGeometrySec[plts_name] = True
                if plts_name in self.U.keys():
                    if break_at_intersection and not np.isclose(self.U[plts_name], 0., 1e-6): 
                        break


            # find last point of the secondary trajectory
            # if (RV_new[0, 0] > 2.5) and (RV_new[0, 1] < 1.5):  
            if (RV_new[0, 0] > r_aim[0] - 0.15) and (RV_new[0, 1] < 1.5):
            # RV_new[0, 4] > 0.0 # Vx > 0
                # intersection with the stop plane:
                r_intersect = gf.plane_segment_intersect(stop_plane_n, r_aim,
                                                      RV_old[0, :3], RV_new[0, :3])
                if r_intersect is not None: 
                    
                    #cut trajectory
                    RV_new[0, :3] = r_intersect
                    RV = np.vstack((RV, RV_new))
                    
                    # check XY plane:
                    if (np.linalg.norm(RV_new[0, :2] - r_aim[:2]) <= eps_xy): #!!!
                        self.IsAimXY = True
                        
                    # check XZ plane:
                    if (np.linalg.norm(RV_new[0, [0, 2]] - r_aim[[0, 2]]) <= eps_z): #!!!
                        self.IsAimZ = True
                    break
                
            # continue trajectory calculation:
            RV_old = RV_new
            t = t + dt
            RV = np.vstack((RV, RV_new))
            tag_column = np.hstack((tag_column, 20))

        else:
            self.print_log("max time exceeded during passing secondary")
            # <krokhalev> debug plot
            plt.figure(348)
            self.plot()
            
        self.RV_sec = RV
        self.tag_sec = tag_column

    def check_B_is_NaN(self, B_local, debug_info=None): 
        if np.isnan(B_local).any():
            self.B_out_of_bounds = True
            self._B_out_of_bounds_debug_info = debug_info #'_pass_to_target'
            # self.print_log('Btor is nan, r = %s' % str(r))
            return True
        else: 
            return False

    def _pass_to_target(self, r_aim, E_interp, B_interp, geom,
                       stop_plane_n=np.array([1., 0., 0.]),
                       eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10,
                       no_intersect=False, no_out_of_bounds=False, 
                       invisible_wall_x=5.5):
        '''
        find secondary trajectory which goes directly to target
        '''
        
        k = self.q / self.m
        
        if True in self.IntersectGeometry.values():
            print('There is intersection at primary trajectory')
            return
        if len(self.Fan) == 0:
            print('NO secondary trajectories')
            return
        
        # reset flags in order to let the algorithm work properly
        self.IsAimXY = False
        self.IsAimZ = False
        
        # reset intersection flags for secondaries
        for key in self.IntersectGeometrySec.keys():
            self.IntersectGeometrySec[key] = False
            
        twisted_fan = self.check_twist(r_aim)
        
        sec_tr_index, linear_param = self.find_fork_for_aim(r_aim, stop_plane_n)
        
        if linear_param is None:
            return
        
        # use runge-kutta to make small step along prim traj
        
        # set start point
        RV_old = np.array([self.Fan[sec_tr_index][0]])
        
        # read fields
        E_local = return_E(RV_old[0, 0:3], E_interp, self.U, geom)
        B_local = B_interp(RV_old[0, 0:3])
        if self.check_B_is_NaN(B_local, (RV_old[0, 0:3], '_pass_to_target')):
             raise Exception("somehow B is None in the middle of primary trajectory")
        
        dt = self.dt1*linear_param
        RV_new = runge_kutt(k, RV_old, dt, E_local, B_local)
        
        # pass new secondary trajectory
        self.pass_sec(RV_new, r_aim, E_interp, B_interp, geom,
                      stop_plane_n=stop_plane_n,
                      eps_xy=eps_xy, eps_z=eps_z, 
                      invisible_wall_x=invisible_wall_x)
        
        # check XY flag
        if self.IsAimXY:
            # insert RV_new into primary traj
            # find the index of the point in primary traj closest to RV_new
            ind_prim = argfind_rv(self.RV_prim, RV_old[0])
            i2insert = ind_prim + 1
            self.RV_prim = np.insert(self.RV_prim, i2insert, RV_new, axis=0)
            self.tag_prim = np.insert(self.tag_prim, i2insert, 11, axis=0)
        else:
            ind_prim = argfind_rv(self.RV_prim, RV_old[0])
            i2insert = ind_prim + 1
            print( "difference along vertical", np.linalg.norm(RV_new[0, :2] - r_aim[:2]) )
            print("linear coeff = ", linear_param)
            print("prim index =", ind_prim)
            print("Fan index =", sec_tr_index)
            plt.figure(555)
            self.plot_prim(plt.gca())
            self.plot_sec(plt.gca())
            self.plot_fan(plt.gca())
            plt.figure(556)
            self.plot_prim(plt.gca(), axes='XZ')
            self.plot_sec(plt.gca(), axes='XZ')
            self.plot_fan(plt.gca(), axes='XZ')
            # raise Exception("didn't make it in the first attempt")
            pass #!!!
    
    def check_twist(self, r_aim):
        
        # find which secondaries are higher/lower than r_aim
        # sign = -1 means higher, 1 means lower
        signs = np.array([np.sign(np.cross(RV[-1, :3], r_aim)[-1])
                          for RV in self.Fan])
        are_higher = np.argwhere(signs == -1)
        are_lower = np.argwhere(signs == 1)
        twisted_fan = False  # flag to detect twist of the fan
        
        try:
            if are_higher[-1] > are_lower[0]:
                print('Fan is twisted!')
                twisted_fan = True
            else:
                self.fan_ok = True
        except IndexError:
            self.print_log("can't check if fan is twisted")
            return None
        
        return twisted_fan
    
    def find_fork_for_aim(self, r_aim, stop_plane_n=np.array([1., 0., 0.])):
        
        '''
        find 2 sec trajectories: 1st is higher than aim, 2nd is lower 
        calc linear parameter t: 0..1 : 0 if r_aim is at 1st traj, 1 if at 2nd
        '''
        
        #set basis
        v, h = gf.vert_horz_basis(stop_plane_n)
        t = None
        idx = None
        
        for i, (rrvv1, rrvv2) in enumerate( zip(self.Fan[0:-1], self.Fan[1:]) ):
            # take last dot
            rv1 = rrvv1[-1, 0:3]
            rv2 = rrvv2[-1, 0:3]
            
            # find vertical difference between last dots and aim
            d1 = (rv1 - r_aim).dot(v)
            d2 = (rv2 - r_aim).dot(v)
            
            # check if last dots on different sides of aim
            if d1*d2 <= 0:
                # calc t
                t = abs(d1/(d1 - d2))
                idx = i
                break
        
        if t is None:
            self.print_log("can't find fork, all secondaries higher or lower than aim")
         
        return idx, t #!!!
            
# %% define class for plates
class Plates():
    
    '''
    object containing info on deflecting plates
    '''

    def __init__(self, name, beamline, r=np.array([0., 0., 0.])):
        '''

        Parameters
        ----------
        name : str
            plates name, 'A1', 'B1' etc
        beamline : str
            beamline name, 'prim' or 'sec'
        r : np.array, optional
            initial plates position. The default is np.array([0., 0., 0.]).

        Returns
        -------
        None.

        '''
        self.name = name
        self.beamline = beamline
        self.r = r
        
        self.rotation_mx = gf.identMx()
        self.inv_rotation_mx = gf.identMx()
        
        # gabarits
        self.min_corners = None
        self.max_corners = None
        
        # direct:  v = plate.rotation_mx.dot(v).A1 + plate.r
        # inverse: v = plate.inv_rotation_mx.dot(v - plate.r).A1 # to original/native coordinates 

    def recalc_gabarits(self, E_interp): 
        g = E_interp[self.name].grid
        edges = [ g[:,  0, 0, 0],  g[:,  0, 0, -1], g[:,  0, -1, 0], g[:,  0, -1, -1], 
                  g[:, -1, 0, 0],  g[:, -1, 0, -1], g[:, -1, -1, 0], g[:, -1, -1, -1] ]
        
        edges_ = np.array( [  self.rotation_mx.dot(e) + self.r for e in edges ] )

        self.min_corners = np.min(edges_, axis=0)
        self.max_corners = np.max(edges_, axis=0)

    def set_edges(self, edges):
        '''
        set coordinates of plates edges
        '''
        self.original_edges = copy.deepcopy(edges)
        self.edges = edges

    def rotate(self, angles, beamline_angles, inverse=False):
        '''
        rotate plates on angles around the axis with beamline_angles
        '''
        self.angles = angles
        self.beamline_angles = beamline_angles

        self.r = gf.rotate3(self.r, angles, beamline_angles, inverse=inverse)  # !!!

        for i in range(self.edges.shape[0]):
            self.edges[i, :] = gf.rotate3(self.edges[i, :], angles, beamline_angles, inverse=inverse)
        
        # remember the transformation 
        mx = gf.get_rotate3_mx(angles, beamline_angles, inverse=inverse)
        self.rotation_mx = mx.dot(self.rotation_mx)
        self.inv_rotation_mx = np.matrix(self.rotation_mx).I.A

    def shift(self, r_new):
        '''
        shift all the coordinates to r_new
        '''
        self.r += r_new
        self.edges += r_new

    def check_intersect(self, point1, point2):
        '''
        check intersection with a segment point1 -> point2
        '''
        #c1, c2 = self.min_corners, self.max_corners
        #if np.any(point1 < c1)and np.any(point2 < c1):
        #    return False
        #if np.any(point1 > c2)and np.any(point2 > c2):
        #    return False

        segment_coords = np.array([point1, point2])
        if gf.segm_poly_intersect(self.edges[0][:4], segment_coords) or \
           gf.segm_poly_intersect(self.edges[1][:4], segment_coords):
            return True
        # if plates are flared
        if self.edges.shape[1] > 4:
            # for flared plates the sequence of vertices is
            # [UP1sw, UP1, UP2, UP2sw, UP3, UP4]
            point_ind = [1, 5, 4, 2]
            if gf.segm_poly_intersect(self.edges[0][point_ind], segment_coords) or \
               gf.segm_poly_intersect(self.edges[1][point_ind], segment_coords):
                return True
            # intersection with flared part
            point_ind = [0, 1, 2, 3]
            if gf.segm_poly_intersect(self.edges[0][point_ind], segment_coords) or \
               gf.segm_poly_intersect(self.edges[1][point_ind], segment_coords):
                return True
        return False

    def plot(self, ax, axes='XY', **kwargs): # <reonid: add **kwargs> #!!!
        '''
        plot plates
        '''
        index_X, index_Y = gf.get_index(axes)
        ax.fill(self.edges[0][:, index_X], self.edges[0][:, index_Y], fill=False, hatch='\\', linewidth=2, **kwargs)
        ax.fill(self.edges[1][:, index_X], self.edges[1][:, index_Y], fill=False, hatch='/', linewidth=2, **kwargs)

    def _rect(self, ii, jj): 
        i0, i1 = ii 
        j0, j1 = jj
        return [ self.edges[i0, j0], self.edges[i1, j0], self.edges[i1, j1], self.edges[i0, j1] ]
        
    def front_rect(self): 
        # checked for A3
        if self.edges.shape[1] > 4:  
            return self._rect((0, 1), (1, 2))
        return [ self.edges[0, 0], self.edges[1, 0], self.edges[1, 1], self.edges[0, 1] ]

    def front_basis(self, norm=True): 
        r = self.front_rect()
        v1 = r[0] - r[1]     # up 
        v2 = r[0] - r[-1]    # right  
        if norm: 
            return [v1/np.linalg.norm(v1),  v2/np.linalg.norm(v2)] # , along_vector
        else: 
            return [v1*0.5,  v2*0.5] 

    def return_r_in_original_coordinates(self, geom, r): 
        # return self.inv_rotation_mx.dot(r - self.r)
        
        r_new = r - geom.r_dict[self.name]
        # get angles
        angles = self.angles
        beamline_angles = self.beamline_angles
        # rotate point to the coord system of plates
        r_new = gf.rotate3(r_new, angles, beamline_angles, inverse=True)
        return r_new
    
    def original_plates(self, geom): 
        _pl = copy.deepcopy(self)
        dr = geom.r_dict[self.name]
        angles = self.angles
        beamline_angles = self.beamline_angles
        
        _pl.shift(-dr)
        _pl.rotate(angles, beamline_angles, inverse=True)
        return _pl        
    
    def contains_point(self, geom, r):
        # !!! not for flare plates 
        if self.edges.shape[1] > 4:  
            raise Exception('for flared plates not implemented yet')
        _r = self.return_r_in_original_coordinates(geom, r)

        _pl = self.original_plates(geom)
               
        xmin = np.min(_pl.edges, axis=(1, 0) )[0]
        ymin = np.min(_pl.edges, axis=(1, 0) )[1]
        zmin = np.min(_pl.edges, axis=(1, 0) )[2]

        xmax = np.max(_pl.edges, axis=(1, 0) )[0]
        ymax = np.max(_pl.edges, axis=(1, 0) )[1]
        zmax = np.max(_pl.edges, axis=(1, 0) )[2]

        _min_r = np.array([xmin, ymin, zmin])
        _max_r = np.array([xmax, ymax, zmax])

        return np.all( _r - _min_r >= 0.0) and np.all( _r - _max_r <= 0.0) 


# %% class for Analyzer
class Analyzer(Plates):
    '''
    Analyzer object
    '''

    def add_slits(self, an_params):
        '''
        add slits and detector to Analyzer
        an_params : list containing [n_slits, slit_dist, slit_w, G, theta,
                                     XD, YD1, YD2]
        n_slits : number of slits
        slit_dist : distance between centers of the slits [m]
        slit_w : slit width (along Y) [m]
        theta : entrance angle to analyzer [deg]
        G : gain function
        XD, YD1, YD2 : geometry parameters [m]
        '''
        # define main parameters of the Analyzer
        self.n_slits, self.slit_dist, self.slit_w, self.G, self.theta, \
            self.XD, self.YD1, self.YD2 = an_params
        self.n_slits = int(self.n_slits)
        # length of the slit
        slit_l = 0.1
        # angles of the slits plane normal
        slit_angles = np.array([self.theta, 0., 0.])
        # coords of center of the central slit
        rs = np.array([0, 0, 0])
        # define slits
        r_slits, slit_plane_n, slits_spot = \
            define_slits(rs, slit_angles, self.n_slits, self.slit_dist,
                         self.slit_w, slit_l)
        # save slits edges
        self.slits_edges = r_slits
        self.slit_plane_n = slit_plane_n
        self.slits_spot = slits_spot
        # define detector
        n_det = self.n_slits
        # set detector angles
        det_angles = np.array([180. - self.theta, 0, 0])
        r_det, det_plane_n, det_spot = \
            define_slits(np.array([self.XD, self.YD1 - self.YD2, 0]),
                         det_angles, n_det, self.slit_dist, self.slit_dist,
                         slit_l)
        # save detector edges
        self.det_edges = r_det
        self.det_plane_n = det_plane_n
        self.det_spot = det_spot
        print('\nAnalyzer with {} slits ok!'.format(self.n_slits))
        print('G = {}'.format(self.G))

    def get_params(self):
        '''
        return analyzer parameters
        [n_slits, slit_dist, slit_w, G, theta, XD, YD1, YD2]
        '''
        print('n_slits = {}\nslit_dist = {}\nslit_width = {}'
              .format(self.n_slits, self.slit_dist, self.slit_w))
        print('G = {}\ntheta = {}\nXD = {}\nYD1 = {}\nYD2 = {}'
              .format(self.G, self.theta, self.XD, self.YD1, self.YD2))
        return(np.array([self.n_slits, self.slit_dist, self.slit_w, self.G,
                         self.theta, self.XD, self.YD1, self.YD2]))

    def rotate(self, angles, beamline_angles, inverse=False):
        '''
        rotate all the coordinates around the axis with beamline_angles
        '''
        super().rotate(angles, beamline_angles, inverse)
        for attr in [self.slits_edges, self.slits_spot,
                     self.det_edges, self.det_spot]:
            if len(attr.shape) < 2:
                attr = gf.rotate3(attr, angles, beamline_angles, inverse)
            else:
                for i in range(attr.shape[0]):
                    attr[i, :] = gf.rotate3(attr[i, :], angles, beamline_angles, inverse)
        # recalculate normal to slit plane:
        self.slit_plane_n = gf.calc_normal(self.slits_edges[0, 0, :],
                                           self.slits_edges[0, 1, :],
                                           self.slits_edges[0, 2, :])
        self.det_plane_n = gf.calc_normal(self.det_edges[0, 0, :],
                                          self.det_edges[0, 1, :],
                                          self.det_edges[0, 2, :])

    def shift(self, r_new):
        '''
        shift all the coordinates to r_new
        '''
        super().shift(r_new)
        for attr in [self.slits_edges, self.slits_spot,
                     self.det_edges, self.det_spot]:
            attr += r_new

    def plot(self, ax, axes='XY', n_slit='all'):
        '''
        plot analyzer
        '''
        # plot plates
        super().plot(ax, axes=axes)
        # choose which slits to plot
        index_X, index_Y = gf.get_index(axes)
        if n_slit == 'all':
            slits = range(self.slits_edges.shape[0])
        else:
            slits = [n_slit]
        # set color cycler
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        colors = colors[:len(slits)]
        colors = cycle(colors)
        # plot slits and detector
        for edges, spot in zip([self.slits_edges, self.det_edges],
                               [self.slits_spot, self.det_spot]):
            for i in slits:
                c = next(colors)
                # plot center
                ax.plot(edges[i, 0, index_X], edges[i, 0, index_Y],
                        '*', color=c)
                # plot edges
                ax.fill(edges[i, 1:, index_X], edges[i, 1:, index_Y],
                        fill=False)
            # plot spot
            ax.fill(spot[:, index_X], spot[:, index_Y], fill=False)
            ax.fill(spot[:, index_X], spot[:, index_Y], fill=False)

def createplates(name, beamline, an_params=None): 
    if name == 'an':
        # create new Analyzer object
        plts = Analyzer(name, beamline)
        plts.add_slits(an_params)
    else:
        plts = Plates(name, beamline)
    return plts

# %%
def define_slits(r0, slit_angles, n_slits, slit_dist, slit_w, slit_l):
    '''
    calculate coordinates of slits edges with central slit at r0
    '''
    n_slits = int(n_slits)
    # calculate slits coordinates:
    r_slits = np.zeros([n_slits, 5, 3])
    for i_slit in range(n_slits):
        # calculate coords of slit center:
        y0 = (n_slits//2 - i_slit)*slit_dist
        r_slits[i_slit, 0, :] = [0., y0, 0.]
        # calculate slit edges:
        r_slits[i_slit, 1, :] = [0., y0 + slit_w/2, slit_l/2]
        r_slits[i_slit, 2, :] = [0., y0 - slit_w/2, slit_l/2]
        r_slits[i_slit, 3, :] = [0., y0 - slit_w/2, -slit_l/2]
        r_slits[i_slit, 4, :] = [0., y0 + slit_w/2, -slit_l/2]
        # rotate and shift to slit position:
        for j in range(5):
            r_slits[i_slit, j, :] = gf.rotate3(r_slits[i_slit, j, :],
                                            slit_angles, slit_angles)
            r_slits[i_slit, j, :] += r0

    # calculate normal to slit plane:
    slit_plane_n = gf.calc_normal(r_slits[0, 0, :], r_slits[0, 1, :],
                               r_slits[0, 2, :])

    # create polygon, which contains all slits (slits spot):
    slits_spot = 1.5*np.vstack([r_slits[0, [1, 4], :] - r0,
                                r_slits[-1, [3, 2], :] - r0]) + r0

    return r_slits, slit_plane_n, slits_spot

# %% define class for geometry
class Geometry():
    '''
    object containing geometry points
    '''

    def __init__(self):
        # dictionary for Plates objects:
        self.plates_dict = dict()
        # dictionary for positions of all objects:
        self.r_dict = dict()
        # arrays for primary and secondary beamline angles:
        self.angles_dict = dict()
        # determining chamber entrance and exit:
        self.chamb_ent = []
        self.chamb_ext = []
        # dictionary for poloidal field coils:
        self.pf_coils = dict()
        # Tor Field coil contour:
        self.coil = np.array([])
        # vacuum vessel contour:
        self.camera = np.array([])
        # separatrix contour:
        self.sep = np.array([])
        # inner and outer first wall contours:
        self.in_fw = np.array([])
        self.out_fw = np.array([])
        # plasma geometry
        self.R = 0
        self.r_plasma = 0
        self.elon = 0

    def check_chamb_intersect(self, beamline, point1, point2):
        '''
        check intersection between segment 1->2 and chamber
        '''
        intersect_flag = False
        # do not check intersection when particle is far from chamber
        if (point1[0] > 2.5 and point2[1] > 1.5) or \
           (point1[0] < 2.0 and point2[1] < 0.8):
               return intersect_flag
        if beamline == 'prim':
            # check intersection with chamber entrance
            # if len(self.chamb_ent) == 0: return False
            for i in np.arange(0, len(self.chamb_ent), 2):
                intersect_flag = intersect_flag or \
                    gf.is_intersect(point1[0:2], point2[0:2],
                                   self.chamb_ent[i], self.chamb_ent[i+1])
        elif beamline == 'sec':
            # check intersection with chamber exit
            # if len(self.chamb_ext) == 0: return False
            for i in np.arange(0, len(self.chamb_ext), 2):
                intersect_flag = intersect_flag or \
                    gf.is_intersect(point1[0:2], point2[0:2],
                                   self.chamb_ext[i], self.chamb_ext[i+1])
        return intersect_flag

    def check_fw_intersect(self, point1, point2):
        '''
        check intersection between segment 1->2 and outer first wall
        '''
        intersect_flag = False
        # do not check intersection when particle is far from first wall
        if (point1[1] > -0.35 and point2[1] > -0.35):
            return intersect_flag
        # check intersection with first wall
        for i in np.arange(4, len(self.out_fw), 2):
            intersect_flag = intersect_flag or \
                gf.is_intersect(point1[0:2], point2[0:2],
                             self.out_fw[i], self.out_fw[i+1])
        return intersect_flag

    def check_plates_intersect(self, point1, point2):
        '''
        check intersection between segment 1->2 and plates
        '''
        # do not check intersection when particle is outside beamlines
        if point2[0] < self.r_dict['aim'][0]-0.05 and point1[1] < self.r_dict['port_in'][1]:
            return False, 'none'
        for key in self.plates_dict.keys():
            # check if a point in inside the beamline
            if (key in ['A1', 'B1', 'A2', 'B2']               and point1[1] > self.r_dict['port_in'][1]) or \
               (key in ['A3', 'A3d', 'B3', 'A4', 'A4d', 'B4'] and point2[0] > self.r_dict['aim'][0]-0.05):
                # check intersection
                if self.plates_dict[key].check_intersect(point1, point2):
                    return True, key
            else:
                continue
        return False, 'none'

    def add_coords(self, name, ref_point, dist, angles):
        '''
        add new element 'name' to r_dict
        '''
        # unpack ref_point
        if type(ref_point) == str:
            r0 = self.r_dict[ref_point]
        else:
            r0 = ref_point
        # unpack angles
        alpha, beta = angles[0:2]
        # coordinates of the center of the object
        r = r0 + gf.calc_vector(dist, alpha, beta)
        self.r_dict[name] = r

    def plot(self, ax, axes='XY', plot_sep=True, plot_aim=True,
             plot_analyzer=True):
        '''
        plot all geometry objects
        '''
        # plot camera and separatrix in XY plane
        if axes == 'XY':
            # plot toroidal coil
            ax.plot(self.coil[:, 0], self.coil[:, 1], '--', color='k')
            ax.plot(self.coil[:, 2], self.coil[:, 3], '--', color='k')
            # plot tokamak camera
            ax.plot(self.camera[:, 0], self.camera[:, 1],
                    color='tab:blue')
            # plot first wall
            ax.plot(self.in_fw[:, 0], self.in_fw[:, 1], color='k')
            ax.plot(self.out_fw[:, 0], self.out_fw[:, 1], color='k')
            # plot separatrix
            if plot_sep:
                ax.plot(self.sep[:, 0] + self.R, self.sep[:, 1],
                        markersize=2, color='b')  # 'tab:orange')
            # plot PF coils
            for coil in self.pf_coils.keys():
                xc = self.pf_coils[coil][0]
                yc = self.pf_coils[coil][1]
                dx = self.pf_coils[coil][2]
                dy = self.pf_coils[coil][3]
                ax.add_patch(Rectangle((xc-dx/2, yc-dy/2), dx, dy,
                                       linewidth=1, edgecolor='tab:gray',
                                       facecolor='tab:gray'))

        index_X, index_Y = gf.get_index(axes)
        # plot plates
        for name in self.plates_dict.keys():
            if name == 'an' and not plot_analyzer:
                continue  # do not plot Analyzer
            self.plates_dict[name].plot(ax, axes=axes)
        if plot_aim:
            # plot aim dot
            ax.plot(self.r_dict['aim'][index_X], self.r_dict['aim'][index_Y],
                    '*', color='b')
            # plot the center of the central slit
            ax.plot(self.r_dict['slit'][index_X], self.r_dict['slit'][index_Y],
                    '*', color='g')

# %%
def add_diafragm(geom, plts_name, diaf_name, diaf_width=0.1):
    '''
    add new plates object which works as a diafragm
    '''
    # create new object in plates dictionary as a copy of existing plates
    geom.plates_dict[diaf_name] = copy.deepcopy(geom.plates_dict[plts_name])
    angles = geom.plates_dict[diaf_name].angles
    beamline_angles = geom.plates_dict[diaf_name].beamline_angles
    r0 = geom.r_dict[plts_name]
    for i in [0, 1]:  # index for upper/lower plate
        for j in [0, 1]:
            # rotate and shift edge to initial coord system
            coords = gf.rotate3(geom.plates_dict[diaf_name].edges[i][j] -
                             r0, angles, beamline_angles, inverse=True)
            # shift up for upper plate and down for lower
            coords += [0, diaf_width*(1-2*i), 0]
            geom.plates_dict[diaf_name].edges[i][3-j] = \
                gf.rotate3(coords, angles, beamline_angles, inverse=False) + r0

# %%
def optimize_B2(tr, geom, UB2, dUB2, E, B, dt, stop_plane_n, target='aim',
                optimize=True, eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10):
    '''
    get voltages on B2 plates and choose secondary trajectory
    which goes into target
    '''
    # set up target
    print('Target: ' + target)
    r_aim = geom.r_dict[target]
    attempts_opt = 0
    attempts_fan = 0
    while True:
        tr.U['B2'], tr.dt1, tr.dt2 = UB2, dt, dt
        # pass fan of secondaries
        tr.pass_fan(r_aim, E, B, geom, stop_plane_n=stop_plane_n,
                    eps_xy=eps_xy, eps_z=eps_z,
                    no_intersect=True, no_out_of_bounds=True, 
                    invisible_wall_x=geom.r_dict[target][0]+0.1)
        # pass trajectory to the target
        tr._pass_to_target(r_aim, E, B, geom, stop_plane_n=stop_plane_n, #!!!
                          eps_xy=eps_xy, eps_z=eps_z, dt_min=dt_min,
                          no_intersect=True, no_out_of_bounds=True, 
                          invisible_wall_x=geom.r_dict[target][0]+0.1)
        print('IsAimXY = ', tr.IsAimXY)
        print('IsAimZ = ', tr.IsAimZ)
        if True in tr.IntersectGeometry.values():
            break
        if not tr.fan_ok:
            attempts_fan += 1
        if attempts_fan > 3 or len(tr.Fan) == 0:
            print('Fan of secondaries is not ok')
            break

        if optimize:
            # change UB2 value proportional to dz
            if not tr.IsAimZ:
                dz = r_aim[2]-tr.RV_sec[-1, 2]
                print('UB2 OLD = {:.2f}, z_aim - z = {:.4f} m'
                      .format(UB2, dz))
                
                UB2_old = UB2 
                UB2 = UB2 - dUB2*dz
                if np.isnan(UB2): 
                    tr.print_log("dUB2 = %f" % dUB2)
                    tr.print_log("dz = %f" % dz)
                    tr.print_log("UB2_old = %f" % UB2_old)
                    
                print('UB2 NEW = {:.2f}'.format(UB2))
                attempts_opt += 1
            else:
                break
            # check if there is a loop while finding secondary to aim
            if attempts_opt > 20:
                print('too many attempts B2!')
                break
        else:
            print('B2 was not optimized')
            break
    return tr

#%% new prim beamline optimization

def single_shot(tr, geom, E, B, dt, stop_plane_n, target='aim', eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10): 
    
        r_aim = geom.r_dict[target]
        
        tr.dt1, tr.dt2 = dt, dt
        # pass fan of secondaries
        #!!! no_intersect=True, no_out_of_bounds=True
        tr.pass_fan       (r_aim, E, B, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z,                no_intersect=False, no_out_of_bounds=False, invisible_wall_x=geom.r_dict[target][0]+0.1)
        # pass trajectory to the target
        tr._pass_to_target(r_aim, E, B, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z, dt_min=dt_min, no_intersect=False, no_out_of_bounds=False, invisible_wall_x=geom.r_dict[target][0]+0.1)
        
        if True in tr.IntersectGeometry.values():
            return None
        
        #if tr.IsAimXY
        dz =  r_aim[2]-tr.RV_sec[-1, 2] # ???
        return dz



def double_shot(tr, geom, UB2_pair, E, B, dt, stop_plane_n, target='aim', eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10): 
    
        _u1, _u2 = UB2_pair 
        tr.U['B2'] = _u1
        dz1 = single_shot(tr, geom,  E, B, dt, stop_plane_n, target='aim', eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10)
        tr.U['B2'] = _u2
        dz2 = single_shot(tr, geom, E, B, dt, stop_plane_n, target='aim', eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10)
        
        if dz1 is None:
            return False, None
        if dz2 is None:
            return False, None
        
        if abs(dz1 -  dz2) < 0.02: 
            if abs(  min(dz1, dz2)  ) < 0.05: 
                return True, _u1
            else: 
                return False, None

        ok = abs(  min(dz1, dz2)  ) < 0.2
                    
        #u_prognosis = _u1 + (_u2 - _u1)*dU_dz(u1, u2, dz1, dz2)
        u_prognosis = _u1 - dz1/(dz2 - dz1)*(_u2 - _u1)
        return ok, u_prognosis
    
def _optimize_B2(tr, geom, UB2, dUB2, E, B, dt, stop_plane_n, target='aim',
                optimize=True, eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10):
    '''
    get voltages on B2 plates and choose secondary trajectory
    which goes into target
    '''
    # set up target
    print('Target: ' + target)
    # r_aim = geom.r_dict[target]


    ok, u_progn = double_shot(tr, geom, (UB2 - 4.0, UB2 + 4.0), E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
    
    if ok: 
        tr.U['B2'] = u_progn
        single_shot(tr, geom, E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
        return tr
    else: 
        if u_progn is not None: 
            ok, u_progn = double_shot(tr, geom, (u_progn - 4.0, u_progn + 4.0), E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
            tr.U['B2'] = u_progn if u_progn is not None else UB2
            single_shot(tr, geom, E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
            return tr
        else: 
            # Brute force
            print('      !!! ----------------- u_progn is None --------------- !!!  ')
            u1, u2 = UB2 - 14.0, UB2 + 14.0  
            ok1, u_progn1 = double_shot(tr, geom, (u1 - 4.0, u1 + 4.0), E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
            ok2, u_progn2 = double_shot(tr, geom, (u2 - 4.0, u2 + 4.0), E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 

            tr.U['B2'] = u_progn1 if u_progn1 is not None else UB2
            dz1 = single_shot(tr, geom, E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
            
            tr.U['B2'] = u_progn2 if u_progn is not None else UB2
            dz2 = single_shot(tr, geom, E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
            
            if dz1 is None:
                return tr
            
            if (dz2 is None) or np.abs(dz1) < np.abs(dz2): 
                tr.U['B2'] = u_progn1 if u_progn1 is not None else UB2
                dz1 = single_shot(tr, geom, E, B, dt, stop_plane_n, target, eps_xy, eps_z, dt_min) 
            return tr

    '''
    #Old parts
    attempts_opt, attempts_fan = 0, 0
    dz, prev_dz = None, None
    ustep, prev_ustep = None, None

    while True:

        tr.U['B2'] = UB2
        tr.dt1, tr.dt2 = dt, dt
        # pass fan of secondaries
        tr.pass_fan       (r_aim, E, B, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z,                no_intersect=True, no_out_of_bounds=True, invisible_wall_x=geom.r_dict[target][0]+0.1)
        # pass trajectory to the target
        tr._pass_to_target(r_aim, E, B, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z, dt_min=dt_min, no_intersect=True, no_out_of_bounds=True, invisible_wall_x=geom.r_dict[target][0]+0.1)
        
        print('IsAimXY = ', tr.IsAimXY); print('IsAimZ = ', tr.IsAimZ)
        if True in tr.IntersectGeometry.values():
            break
        if not tr.fan_ok:
            attempts_fan += 1
        if attempts_fan > 3 or len(tr.Fan) == 0:
            print('Fan of secondaries is not ok')
            break

        # change UB2 value proportional to dz
        if not tr.IsAimZ:
            prev_dz, dz = dz, r_aim[2]-tr.RV_sec[-1, 2]
            if (prev_dz is not None)and(prev_ustep is not None): 
                print('before    dz=', dz, '   ustep=', ustep, '   prev_ustep=', prev_ustep, '    prev_dz=', prev_dz)
                prev_ustep, ustep = ustep, prev_ustep*dz/(dz - prev_dz)
                print('after    dz=', dz, '   ustep=', ustep, '   prev_ustep=', prev_ustep, '    prev_dz=', prev_dz)
            else:     
                print('before    dz=', dz, '   ustep=', ustep, '   prev_ustep=', prev_ustep, '    prev_dz=', prev_dz)
                prev_ustep, ustep = ustep, - dUB2*dz
                print('after    dz=', dz, '   ustep=', ustep, '   prev_ustep=', prev_ustep, '    prev_dz=', prev_dz)
             
            UB2 = UB2 + ustep
            attempts_opt += 1

        else:
            break
        
        if attempts_opt > 20:  # check if there is a loop while finding secondary to aim
            print('too many attempts B2!')
            break
    return tr
    '''

# %%
def calc_zones(tr, dt, E, B, geom, slits=[2], timestep_divider=5,
               stop_plane_n=np.array([1., 0., 0.]), eps_xy=1e-3, eps_z=1,
               dt_min=1e-11, no_intersect=True, no_out_of_bounds=True):
    '''
    calculate ionization zones
    '''
    # find the number of slits
    n_slits = geom.plates_dict['an'].slits_edges.shape[0]
    tr.add_slits(n_slits)
    # set target at the central slit
    r_aim = geom.plates_dict['an'].slits_edges[n_slits//2, 0, :]

    # create slits polygon
    slit_plane_n = geom.plates_dict['an'].slit_plane_n
    slits_spot = geom.plates_dict['an'].slits_spot
    ax_index = np.argmax(slit_plane_n)
    slits_spot_flat = np.delete(slits_spot, ax_index, 1)
    slits_spot_poly = path.Path(slits_spot_flat)

    # find trajectories which go to upper and lower slit edge
    # find index of primary trajectory point where secondary starts
    index = np.nanargmin(np.linalg.norm(tr.RV_prim[:, :3] -
                                        tr.RV_sec[0, :3], axis=1))
    # set up the index range
    sec_ind = range(index-2, index+2)
    print('\nStarting precise fan calculation')
    k = tr.q / tr.m
    # divide the timestep
    tr.dt1 = dt/timestep_divider
    tr.dt2 = dt
    # number of steps during new fan calculation
    n_steps = timestep_divider * (len(sec_ind))
    # list for new trajectories
    fan_list = []
    # take the point to start fan calculation
    # RV_old = tr.Fan[sec_ind[0]-1][0]
    RV_old = tr.RV_prim[sec_ind[0]]
    RV_old = np.array([RV_old])
    RV_new = RV_old

    i_steps = 0
    while i_steps <= n_steps:
        # pass new secondary trajectory
        tr.pass_sec(RV_new, r_aim, E, B, geom,
                    stop_plane_n=slit_plane_n,
                    tmax=9e-5, eps_xy=1e-3, eps_z=1)
        
        # make a step on primary trajectory
        r = RV_old[0, :3]
        
        #fields
        E_local = np.array([0., 0., 0.])
        B_local = B(r)
        if np.isnan(B_local).any():
            # print('Btor is nan, r = %s' % str(r))
            break
                
        #runge-kutta step
        RV_new = runge_kutt(k, RV_old, tr.dt1, E_local, B_local)
        RV_old = RV_new
        i_steps += 1
        
        # check intersection with slits polygon
        intersect_coords_flat = np.delete(tr.RV_sec[-1, :3], ax_index, 0)
        contains_point = slits_spot_poly.contains_point(intersect_coords_flat)
        if not (True in tr.IntersectGeometrySec.values() or
                tr.B_out_of_bounds) and contains_point:
            fan_list.append(tr.RV_sec)
    tr.Fan = fan_list
    tr.fan_to_slits = fan_list
    print('\nPrecise fan calculated')

    for i_slit in slits:
        # set up upper and lower slit edge
        upper_edge = [geom.plates_dict['an'].slits_edges[i_slit, 4, :],
                      geom.plates_dict['an'].slits_edges[i_slit, 1, :]]
        lower_edge = [geom.plates_dict['an'].slits_edges[i_slit, 3, :],
                      geom.plates_dict['an'].slits_edges[i_slit, 2, :]]
        zones_list = []  # list for ion zones coordinates
        rv_list = []  # list for RV arrays of secondaries
        for edge in [upper_edge, lower_edge]:
            # find intersection of fan and slit edge
            for i_tr in range(len(tr.Fan) - 1):
                p1, p2 = tr.Fan[i_tr][-1, :3], tr.Fan[i_tr+1][-1, :3]
                # check intersection between fan segment and slit edge
                if gf.is_intersect(p1, p2, edge[0], edge[1]):
                    r_intersect = gf.segm_intersect(p1, p2, edge[0], edge[1])
                    print('\n intersection with slit ' + str(i_slit))
                    print(r_intersect)
                    tr.dt1 = dt/timestep_divider
                    tr.dt2 = dt
                    tr.pass_to_target(r_intersect, E, B, geom,
                                      eps_xy=eps_xy, eps_z=eps_z,
                                      dt_min=dt_min,
                                      stop_plane_n=slit_plane_n,
                                      no_intersect=no_intersect,
                                      no_out_of_bounds=no_out_of_bounds)
                    zones_list.append(tr.RV_sec[-1, :3])
                    rv_list.append(tr.RV_sec)
                    print('ok!')
                    break
        tr.ion_zones[i_slit] = np.array(zones_list)
        tr.RV_sec_toslits[i_slit] = rv_list
    return tr


# %%
def pass_to_slits(tr, dt, E, B, geom, target='slit', timestep_divider=10,
                   slits=range(5), no_intersect=True, no_out_of_bounds=True):
    '''
    pass trajectories to slits and save secondaries which get into slits
    '''
    tr.dt1 = dt/4
    tr.dt2 = dt
    k = tr.q / tr.m
    # find the number of slits
    n_slits = geom.plates_dict['an'].slits_edges.shape[0]
    tr.add_slits(n_slits)
    # find slits position
    if target == 'slit':
        r_slits = geom.plates_dict['an'].slits_edges
        rs = geom.r_dict['slit']
        slit_plane_n = geom.plates_dict['an'].slit_plane_n
        slits_spot = geom.plates_dict['an'].slits_spot
    elif target == 'det':
        r_slits = geom.plates_dict['an'].det_edges
        rs = geom.r_dict['det']
        slit_plane_n = geom.plates_dict['an'].det_plane_n
        slits_spot = geom.plates_dict['an'].det_spot

    # create slits polygon
    ax_index = np.argmax(slit_plane_n)
    slits_spot_flat = np.delete(slits_spot, ax_index, 1)
    slits_spot_poly = path.Path(slits_spot_flat)

    # find index of primary trajectory point where secondary starts
    index = np.nanargmin(np.linalg.norm(tr.RV_prim[:, :3] -
                                        tr.RV_sec[0, :3], axis=1))
    sec_ind = range(index-2, index+2)

    print('\nStarting precise fan calculation')
    # divide the timestep
    tr.dt1 = dt/timestep_divider
    tr.dt2 = dt
    # number of steps during new fan calculation
    n_steps = timestep_divider * (len(sec_ind))
    # list for new trajectories
    fan_list = []
    # take the point to start fan calculation
    # RV_old = tr.Fan[sec_ind[0]-1][0]
    RV_old = tr.RV_prim[sec_ind[0]]
    RV_old = np.array([RV_old])
    RV_new = RV_old

    i_steps = 0
    inside_slits_poly = False
    while i_steps <= n_steps:
        # pass new secondary trajectory
        tr.pass_sec(RV_new, rs, E, B, geom,
                    stop_plane_n=slit_plane_n, tmax=9e-5,
                    eps_xy=1e-3, eps_z=1)
        
        # make a step on primary trajectory
        r = RV_old[0, :3]
        
        #fields
        E_local = np.array([0., 0., 0.])
        B_local = B(r)
        if np.isnan(B_local).any():
            print('Btor is nan, r = %s' % str(r))
            break
        
        #runge-kutta step
        RV_new = runge_kutt(k, RV_old, tr.dt1, E_local, B_local)
        RV_old = RV_new
        i_steps += 1
        
        # ???
        intersect_coords_flat = np.delete(tr.RV_sec[-1, :3], ax_index, 0)
        contains_point = slits_spot_poly.contains_point(intersect_coords_flat)
        
        # save result
        if not (True in tr.IntersectGeometrySec.values() or
                tr.B_out_of_bounds) and contains_point:
            inside_slits_poly = True
            fan_list.append(tr.RV_sec)
        if not contains_point and inside_slits_poly:
            break
    print('\nPrecise fan calculated')

    # choose secondaries which get into slits
    # start slit cycle
    for i_slit in slits:
        print('\nslit = {}'.format(i_slit+1))
        print('center of the slit = ', r_slits[i_slit, 0, :], '\n')
        
        # create slit polygon
        slit_flat = np.delete(r_slits[i_slit, 1:, :], ax_index, 1)
        slit_poly = path.Path(slit_flat)
        zones_list = []  # list for ion zones coordinates
        rv_list = []  # list for RV arrays of secondaries
        
        for fan_tr in fan_list:
            # get last coordinates of the secondary trajectory
            intersect_coords_flat = np.delete(fan_tr[-1, :3], ax_index, 0)
            if slit_poly.contains_point(intersect_coords_flat):
                print('slit {} ok!\n'.format(i_slit+1))
                rv_list.append(fan_tr)
                zones_list.append(fan_tr[0, :3])

        tr.RV_sec_toslits[i_slit] = rv_list
        tr.ion_zones[i_slit] = np.array(zones_list)
    tr.fan_to_slits = fan_list

    return tr

# %%
def save_traj_list(traj_list, Btor, Ipl, r_aim, dirname='output'):
    '''
    save list of Traj objects to *.pkl file
    '''

    if len(traj_list) == 0:
        print('traj_list empty! nothing to save')
        return

    Ebeam_list = []
    UA2_list = []

    for traj in traj_list:
        Ebeam_list.append(traj.Ebeam)
        UA2_list.append(traj.U['A2'])

    dirname = dirname + '/' + 'B{}_I{}'.format(int(Btor), int(Ipl))

    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname, 0o700)
            print("Directory ", dirname, " created ")
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    fname = dirname + '/' + 'E{}-{}'.format(int(min(Ebeam_list)),
                                            int(max(Ebeam_list))) + \
        '_UA2{}-{}'.format(int(min(UA2_list)), int(max(UA2_list))) + \
        '_alpha{:.1f}_beta{:.1f}'.format(traj.alpha, traj.beta) +\
        '_x{}y{}z{}.pkl'.format(int(r_aim[0]*100), int(r_aim[1]*100),
                                int(r_aim[2]*100))

    with open(fname, 'wb') as f:
        pc.dump(traj_list, f, -1)

    print('\nSAVED LIST: \n' + fname)

# %%
def read_traj_list(fname, dirname='output'):
    '''
    import list of Traj objects from *.pkl file
    '''
    with open(dirname + '/' + fname, 'rb') as f:
        traj_list = pc.load(f)
    return traj_list

# %%
def save_traj2dat(traj_list, save_fan=False, dirname='output/',
                  fmt='%.2f', delimiter=' '):
    '''
    save list of trajectories to *.dat files for CATIA plot
    '''
    for tr in traj_list:
        # save primary
        fname = dirname + 'E{:.0f}_U{:.0f}_prim.dat'.format(tr.Ebeam, tr.U['A2'])
        np.savetxt(fname, tr.RV_prim[:, 0:3]*1000,
                   fmt=fmt, delimiter=delimiter)  # [mm]
        # save secondary
        fname = dirname + 'E{:.0f}_U{:.0f}_sec.dat'.format(tr.Ebeam, tr.U['A2'])
        np.savetxt(fname, tr.RV_sec[:, 0:3]*1000,
                   fmt=fmt, delimiter=delimiter)

#%% sec beamline optimization

def optimize(tr, optimizer, E, B, geom, RV0, tmax, eps_xy, eps_z):
    while True:
        # print('\n passing secondary trajectory')
        tr.pass_sec(RV0, optimizer.r_target, E, B, geom, tmax=tmax, 
                    eps_xy=eps_xy, eps_z=eps_z, stop_plane_n=optimizer.stop_plane_n)
        if not optimizer.change_voltages(tr): 
            return optimizer.voltages_list
        # print('\n changing voltages')

def _optimize(tr, optimizer, E, B, geom, RV0, tmax, eps_xy, eps_z, plot_all=False):
    uu = np.linspace(-40.0, 40.0, 41)
    dd = np.zeros_like(uu)
    # dd_v, dd_h, dd_x = np.zeros_like(uu), np.zeros_like(uu), np.zeros_like(uu)
    for i, u in enumerate( uu ): 
        optimizer.set_voltage(tr, u)
        tr._pass_sec(RV0, optimizer.r_target, E, B, geom, tmax=tmax, 
                eps_xy=eps_xy, eps_z=eps_z, stop_plane_n=optimizer.stop_plane_n, 
                break_at_intersection=True)
        dd[i] = optimizer.calc_delta(tr.RV_sec[-1])
 
                                               #dd_v[i], dd_h[i] = optimizer.calc_delta_vector(last_r)
        # dd_x[i], dd_v[i], dd_h[i] = last_r
        if plot_all:
            plt.figure(225)
            tr.plot_prim(plt.gca(), 'XY')
            tr.plot_sec(plt.gca(), 'XY')
    
            plt.figure(229)
            tr.plot_prim(plt.gca(), 'XZ')
            tr.plot_sec(plt.gca(), 'XZ')
        # print(dd[i])
    
    # plt.figure(223)
    # plt.plot(uu, dd)

            plt.figure(300)
            plt.plot(uu, dd)
    # plt.figure(301)
    # plt.plot(uu, dd_h)
    # plt.figure(302)
    # plt.plot(dd_h, dd_v)
    # plt.plot(dd_x, dd_h)

    
    i, iplus1, t = find_fork(dd, threshold=optimizer.aim_zone_size)  
    if i is None: 
        return []
    else: 
        U =  uu[i] + (uu[iplus1]-uu[i])*t 
        optimizer.set_voltage(tr, U)
        tr._pass_sec(RV0, optimizer.r_target, E, B, geom, tmax=tmax, 
                eps_xy=eps_xy, eps_z=eps_z, stop_plane_n=optimizer.stop_plane_n, 
                break_at_intersection=True)

        return [U]
        

        # print('\n changing voltages')

def optimize_sec_fine(E, B, geom, traj_list, optimization_mode="center", 
                 target='slit', max_voltages=[40., 40., 40.], eps_xy=1e-3, 
                 eps_z=1e-3, tmax=9e-5, U0 = [0., 0., 0.], dU = [7, 10, 5.]):
    
    #create list for results
    traj_list_passed = []
    
    #create 3 optimizers
    opt_A3 = optimizers.Optimizer('zone',            dU[0], max_voltages[0], 'A4', geom, 'A3', aim_rough=0.2)
    print('aim_zone_size', opt_A3.aim_zone_size)
    opt_B3 = optimizers.Optimizer(optimization_mode, dU[1], max_voltages[1], target, geom, 'B3')
    opt_A4 = optimizers.Optimizer(optimization_mode, dU[2], max_voltages[2], target, geom, 'A4')
    opt_B3_A4 = optimizers.Multiple_optimizer([opt_B3, opt_A4])
    
    for tr in traj_list:
        
        print('\nEb = {}, UA2 = {}'.format(tr.Ebeam, tr.U['A2']))
        print('Target: ' + target)
        
        #set start point and voltages
        RV0 = np.array([tr.RV_sec[0]])
        tr.U['A3'], tr.U['B3'], tr.U['A4'] = U0
        
        #reset optimizers
        A3_voltages, B3_voltages, A4_voltages = [], [], []
        opt_A3.voltages_list = []
        opt_A3.count = 0
        opt_A3.U_step= None
        for optimizer in opt_B3_A4.optimizers:
            optimizer.voltages_list = []
            optimizer.count = 0
            optimizer.U_step= None
        
        # -------------- optimize plates one by one --------------------
        A3_voltages = _optimize(tr, opt_A3, E, B, geom, RV0, tmax, eps_xy, eps_z, plot_all=True)
        B3_voltages = _optimize(tr, opt_B3, E, B, geom, RV0, tmax, eps_xy, eps_z)
        A4_voltages = _optimize(tr, opt_A4, E, B, geom, RV0, tmax, eps_xy, eps_z)
        print("A3 voltages: ", A3_voltages)
        print("B3 voltages: ", B3_voltages)
        print("A4 voltages: ", A4_voltages)
        
        #finish optimization and write down results
        if (not A4_voltages) or (not B3_voltages) or (True in tr.IntersectGeometrySec.values()):
            print("\nOptimization failed, trajectory NOT saved\n")
        else: 
            print("\nTrajectory optimized and saved\n")
            traj_list_passed.append(tr)
            
    return traj_list_passed

