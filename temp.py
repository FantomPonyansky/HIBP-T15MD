class Traj_():
    ...
    
        def check_intersect_prim(self, geom, RV_old, RV_new, invisible_wall_x): 
        #check if out of bounds for passing to aim
        if RV_new[0, 0] > invisible_wall_x and RV_new[0, 1] < 1.2:
            return True  # self.print_log('primary hit invisible wall, r = %s' % str(np.round(r, 3)))

        if geom.check_chamb_intersect('prim', RV_old[0, 0:3], RV_new[0, 0:3]):
            self.IntersectGeometry['chamb'] = True
            return True  # self.print_log('Primary intersected chamber entrance')

        plts_flag, plts_name = geom.check_plates_intersect(RV_old[0, 0:3], RV_new[0, 0:3])
        if plts_flag:
            self.IntersectGeometry[plts_name] = True
            return True   # #self.print_log('Primary intersected ' + plts_name + ' plates')

        if geom.check_fw_intersect(RV_old[0, 0:3], RV_new[0, 0:3]): 
            return True  # stop primary trajectory calculation  # self.print_log('Primary intersected first wall')

        return False

    #%%----------------------------------------------------------------------------
    def pass_prim(self, E_interp, B_interp, geom, tmax=1e-5, 
                  invisible_wall_x=5.5):
        '''
        passing primary trajectory from initial point self.RV0
        E_interp : dictionary with E field interpolants
        B_interp : list with B fied interpolants
        geom : Geometry object
        '''
        self.reset_flags("only_intersection_flags")
        t = 0.
        dt = self.dt1
        RV_old = self.RV0  # initial position
        RV = self.RV0  # array to collect all r, V
        k = self.q / self.m
        tag_column = [10]

        while t <= tmax:
            r = RV_old[0, :3]
            E_local = return_E(r, E_interp, self.U, geom)
            B_local = return_B(r, B_interp); 
            if np.isnan(B_local).any(): break # self.print_log('Btor is nan, r = %s' % str(r))
            # *********************************************************************
            RV_new = runge_kutt(k, RV_old, dt, E_local, B_local)
            # *********************************************************************
            RV, tag_column =  np.vstack((RV, RV_new)),  np.hstack((tag_column, 10))
            if self.check_intersect_prim(geom, RV_old, RV_new, invisible_wall_x): break  
            RV_old, t = RV_new, t + dt
            
        self.RV_prim, self.tag_prim = RV, tag_column
        
    #-----------------------------------------------------------------------------
    
    def check_aim(self, RV_new, r_aim, eps_xy=1e-3, eps_z=1e-3): 
        # check XY plane:    
        if (np.linalg.norm(RV_new[0, :2] - r_aim[:2]) <= eps_xy): 
            self.IsAimXY = True   # print('aim XY!')
        # check XZ plane:
        if (np.linalg.norm(RV_new[0, [0, 2]] - r_aim[[0, 2]]) <= eps_z): 
            self.IsAimZ = True    # print('aim Z!')    

def ___No_break___(): 
    pass

    def check_intersect_sec(self, geom, RV_old, RV_new, invisible_wall_x): 
        # #check if out of bounds for passing to aim
        if RV_new[0, 0] > invisible_wall_x:
            return True, None  # self.print_log('secondary hit invisible wall, r = %s' % str(np.round(r, 3)))
            
        if geom.check_chamb_intersect('sec', RV_old[0, 0:3], RV_new[0, 0:3]):
            self.IntersectGeometrySec['chamb'] = True  # print('Secondary intersected chamber exit')
            ___No_break___ # ??? Why

        plts_flag, plts_name = geom.check_plates_intersect(RV_old[0, 0:3], RV_new[0, 0:3])
        if plts_flag:
            self.IntersectGeometrySec[plts_name] = True # self.print_log('Secondary intersected ' + plts_name + ' plates')
            ___No_break___ # ??? Why

        # find last point of the secondary trajectory
        if (RV_new[0, 0] > 2.5) and (RV_new[0, 1] < 1.5): # if Last point is outside
            # intersection with the stop plane:
            r_intersect = line_plane_intersect(stop_plane_n, r_aim, RV_new[0, :3]-RV_old[0, :3], RV_new[0, :3])
            
            # check if r_intersect is between RV_old and RV_new:
            if is_between(RV_old[0, :3], RV_new[0, :3], r_intersect):
                #RV_new[0, :3] = r_intersect
                #RV = np.vstack((RV, RV_new))
                # tag_column = np.hstack((tag_column, ??)) ????????? Не добавлен таг. Забыт или так надо?  
                self.IsAimXY, self.IsAimZ = self.check_aim(RV_new, r_aim, eps_xy, eps_z)  # check XY plane, check XZ plane
                return True, r_intersect
          return False, None  


    #%%-----------------------------------------------------------------------------
    def pass_sec(self, RV0, r_aim, E_interp, B_interp, geom,
                 stop_plane_n=np.array([1., 0., 0.]), tmax=5e-5,
                 eps_xy=1e-3, eps_z=1e-3, invisible_wall_x=5.5):
        '''
        passing secondary trajectory from initial point RV0 to point r_aim
        with accuracy eps
        RV0 : initial position and velocity
        '''
        # print('Passing secondary trajectory')
        self.reset_flags("all")

        t = 0.
        dt = self.dt2
        RV_old = RV0  # initial position
        RV = RV0  # array to collect all [r,V]
        tag_column = [20]
        k = 2*self.q / self.m

        while t <= tmax:  # to witness curls initiate tmax as 1 sec
            r = RV_old[0, :3]
            E_local = return_E(r, E_interp, self.U, geom)
            B_local = return_B(r, B_interp); 
            if np.isnan(B_local).any(): break  # self.print_log('Btor is nan, r = %s' % str(r))
            # ******************************************************
            RV_new = runge_kutt(k, RV_old, dt, E_local, B_local)
            # ****************************************************** 
        
            stop, r_intersect = self.check_intersect_sec(geom, RV_old, RV_new, invisible_wall_x)
            if r_intersect in not None: 
                RV_new[0, :3] = r_intersect
                RV = np.vstack((RV, RV_new))
                # tag_column = np.hstack((tag_column, ??)) ????????? Не добавлен таг. Забыт или так надо?                  
            if stop: 
                break
        
            # continue trajectory calculation:
            RV_old, t = RV_new, t + dt
            RV, tag_column = np.vstack((RV, RV_new)),   np.hstack((tag_column, 20))

        self.RV_sec, self.tag_sec = RV, tag_column

    #%%-----------------------------------------------------------------------------
            
    def pass_fan(self, r_aim, E_interp, B_interp, geom,
                 stop_plane_n=np.array([1., 0., 0.]), eps_xy=1e-3, eps_z=1e-3,
                 no_intersect=False, no_out_of_bounds=False, 
                 invisible_wall_x=5.5):
        '''
        passing fan from initial point self.RV0
        '''
        # ********************************************************* #               
        self.pass_prim(E_interp, B_interp, geom, invisible_wall_x=invisible_wall_x)
        # ********************************************************* #               
                       
        # create a list fro secondary trajectories:
        list_sec = []
        # check intersection of primary trajectory:
        if True in self.IntersectGeometry.values():
            self.Fan = [] # list_sec   # print('Fan list is empty')
            return

        # check eliptical radius of particle:  # 1.5 m - major radius of a torus, elon - size along Y
        mask = np.sqrt((self.RV_prim[:, 0] - geom.R)**2 + (self.RV_prim[:, 1] / geom.elon)**2) <= geom.r_plasma   
        
        
#%%
if not load_traj_from_file:
    # define list of trajectories that hit r_aim
    traj_list_B2 = []
    # initial beam energy range
    Ebeam_range = np.arange(Emin, Emax + dEbeam, dEbeam)  # [keV]
    
    for Ebeam in Ebeam_range:
        shift, geomT15.r_dict['aim_zshift'], dUB2, t1, shot, input_fname, exp_voltages, indexes = ...
        target, UA2_range, UB2_range, UA3_range, UB3_range, optimizeA3B3, eps_xy, eps_z = ...
    
        # UA2 loop
        for i in range(UA2_range.shape[0]):
            UA2 = UA2_range[i]
            if not optimizeB2: UB2 = UB2_range[i]
            if not optimizeA3B3: UA3, UB3 = UA3_range[i], UB3_range[i]
                
            #reset aim point
            shift = np.zeros(3)
            geomT15.r_dict['aim_zshift'] = geomT15.r_dict['aim']
            
            # dict of starting voltages
            U_dict = {'A2': UA2, 'B2': UB2,'A3': UA3, 'B3': UB3, 'A4': UA4, 'an': Ebeam/(2*G)}
            tr = hb.Traj(q, m_ion, Ebeam, r0, geomT15.angles_dict['r0'][0], geomT15.angles_dict['r0'][1], U_dict, dt)
            
            # **************************** optimize B2 voltage. here the trajectories calculated !!!
            tr = hb.optimize_B2(tr, geomT15, UB2, dUB2, E, B, dt, stop_plane_n, target, optimizeB2, eps_xy=eps_xy, eps_z=eps_z)
            if tr is None: continue #if True in tr.IntersectGeometry.values()  if True in tr.IntersectGeometrySec.values(): continue # print('NOT saved, primary intersected geometry')
           # **************************** 
            
 
            # if no intersections, upldate UB2 values
            UB2 = tr.U['B2']

            #calc shift
            shift = calc_shift(tr.RV_sec[-1, 3:6], geomT15.plates_dict['A3'], adaptive_aim)
            if is_zero(shift):  # np.all( np.isclose(shift, np.zeros(3)) ): 
                traj_list_B2.append(tr)
            else:
                geomT15.r_dict['aim_zshift'] = geomT15.r_dict['aim'] + shift
                tr_shifted = hb.optimize_B2(tr, geomT15, UB2, dUB2, E, B, dt, stop_plane_n, target, optimizeB2, eps_xy=eps_xy, eps_z=eps_z)
                traj_list_B2.append(     best(shifted_tr, tr, ...)     )

            else:
                print('NOT saved, sth is wrong')
    
    t2 = time.time()
    print({True: '\n B2 voltage optimized, t = {:.1f} s\n'.format(t2-t1), False: '\n Trajectories to r_aim calculated, t = {:.1f} s\n'.format(t2-t1)}[optimizeB2])

#%%  pure experiment
if not load_traj_from_file:
    # define list of trajectories that hit r_aim
    traj_list_B2 = []
    # initial beam energy range
    Ebeam_range = np.arange(Emin, Emax + dEbeam, dEbeam)  # [keV]
    
    for Ebeam in Ebeam_range:
        shift, geomT15.r_dict['aim_zshift'], dUB2, t1, shot, input_fname, exp_voltages, indexes = ...
        target, UA2_range, UB2_range, UA3_range, UB3_range, optimizeA3B3, eps_xy, eps_z = ...
    
        # UA2 loop
        for i in range(UA2_range.shape[0]):
            UA2 = UA2_range[i]
            UB2 = UB2_range[i]
            UA3, UB3 = UA3_range[i], UB3_range[i]
                
            #reset aim point
            shift = np.zeros(3)
            geomT15.r_dict['aim_zshift'] = geomT15.r_dict['aim']
            
            # dict of starting voltages
            U_dict = {'A2': UA2, 'B2': UB2,'A3': UA3, 'B3': UB3, 'A4': UA4, 'an': Ebeam/(2*G)}
            tr = hb.Traj(q, m_ion, Ebeam, r0, geomT15.angles_dict['r0'][0], geomT15.angles_dict['r0'][1], U_dict, dt)
            
            # **************************** optimize B2 voltage. here the trajectories calculated !!!
            tr = hb.optimize_B2(tr, geomT15, UB2, dUB2, E, B, dt, stop_plane_n, target, False, eps_xy=eps_xy, eps_z=eps_z)
            if tr is None: continue #if True in tr.IntersectGeometry.values()  if True in tr.IntersectGeometrySec.values(): continue # print('NOT saved, primary intersected geometry')
            # **************************** 
            traj_list_B2.append(tr)
    
    t2 = time.time()
    print({True: '\n B2 voltage optimized, t = {:.1f} s\n'.format(t2-t1), False: '\n Trajectories to r_aim calculated, t = {:.1f} s\n'.format(t2-t1)}[optimizeB2])


#%%

        
        self.tag_prim[mask] = 11

        # list of initial points of secondary trajectories:
        RV0_sec = self.RV_prim[(self.tag_prim == 11)]

        for RV02 in RV0_sec:
            RV02 = np.array([RV02])
            # ********************************************************* #
            self.pass_sec(RV02, r_aim, E_interp, B_interp, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z, invisible_wall_x=invisible_wall_x)
            # ********************************************************* #
            if sec_traj_is_ok(self, no_intersect, no_out_of_bounds): # not (     (no_intersect and True in self.IntersectGeometrySec.values()) or (no_out_of_bounds and self.B_out_of_bounds)    ):
                list_sec.append(self.RV_sec)


        self.Fan = list_sec
        
    #%%-----------------------------------------------------------------------------


    def check_fan_twisted(self, r_aim):
        signs = np.array([np.sign(np.cross(RV[-1, :3], r_aim)[-1]) for RV in self.Fan])
        are_higher = np.argwhere(signs == -1)
        are_lower = np.argwhere(signs == 1)
        twisted_fan = False  # flag to detect twist of the fan

        if are_higher.shape[0] == 0:
            n = int(are_lower[are_lower.shape[0]//2]) # print('all secondaries are lower than aim!')
        elif are_lower.shape[0] == 0:
            n = int(are_higher[are_higher.shape[0]//2]) # print('all secondaries are higher than aim!')
        else:
            if are_higher[-1] > are_lower[0]:
                twisted_fan = True  # print('Fan is twisted!')
                n = int(are_lower[-1])
            else:
                n = int(are_higher[-1])  # find the last one which is higher
                self.fan_ok = True
        return n, twisted_fan

    def reset_flags(self, options): 
        if (options == "all") or (options == "intersection_flags_and_aim"): 
            self.IsAimXY = False
            self.IsAimZ = False
            
        if (options == "all"): 
            self.B_out_of_bounds = False
        
        if (options == "only_intersection_flags") or (options == "all") or (options == "intersection_flags_and_aim"): 
            # reset intersection flags for secondaries
            for key in self.IntersectGeometrySec.keys():
                self.IntersectGeometrySec[key] = False
    
    def skip_pass_to_target(self):
        if True in self.IntersectGeometry.values(): 
            return True # print('There is intersection at primary trajectory'); 
        elif len(self.Fan) == 0: 
            return True # print('NO secondary trajectories'); 
        else:
            return False        
    
    def find_the_index_of_the_point_in_primary_traj_closest_to(self, RV_new): 
        # insert RV_new into primary traj
        # find the index of the point in primary traj closest to RV_new
        ind = np.nanargmin(np.linalg.norm(self.RV_prim[:, :3] - RV_new[0, :3], axis=1))
        if is_between(self.RV_prim[ind, :3], self.RV_prim[ind+1, :3], RV_new[0, :3], eps=1e-4):
            i2insert = ind+1
        else:
            i2insert = ind    
        return i2insert 
    
    #%%-----------------------------------------------------------------------------
    def pass_to_target(self, r_aim, E_interp, B_interp, geom,
                       stop_plane_n=np.array([1., 0., 0.]),
                       eps_xy=1e-3, eps_z=1e-3, dt_min=1e-10,
                       no_intersect=False, no_out_of_bounds=False, 
                       invisible_wall_x=5.5):
        '''
        find secondary trajectory which goes directly to target
        '''
        if self.skip_pass_to_target(): return  
        self.reset_flags("intersection_flags_and_aim")

        # find which secondaries are higher/lower than r_aim   # sign = -1 means higher, 1 means lower        
        n, twisted = self.check_fan_twisted(r_aim)  
        RV_old = np.array([self.Fan[n][0]])

        # find secondary, which goes directly into r_aim
        self.dt1 = self.dt1/2.
        while True:
            # make a small step along primary trajectory
            r = RV_old[0, :3]
            B_local = return_B(r, B_interp);   
            if np.isnan(B_local).any(): break
            E_local = np.array([0., 0., 0.])
            # ********************************************************* #
            RV_new = runge_kutt(self.q / self.m, RV_old, self.dt1, E_local, B_local)                                   
            # pass new secondary trajectory
            # ********************************************************* #
            self.pass_sec(RV_new, r_aim, E_interp, B_interp, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z, invisible_wall_x=invisible_wall_x)
            # ********************************************************* #               

            # check XY flag
            if self.IsAimXY:
                # insert RV_new into primary traj
                i2insert = self.find_the_index_of_the_point_in_primary_traj_closest_to(RV_new)
                self.RV_prim, self.tag_prim = np.insert(self.RV_prim, i2insert, RV_new, axis=0),   np.insert(self.tag_prim, i2insert, 11, axis=0)
                break

            # check if the new secondary traj is lower than r_aim
            if (not twisted and np.sign(np.cross(self.RV_sec[-1, :3], r_aim)[-1]) > 0): # if lower, halve the timestep and try once more
                self.dt1 = self.dt1/2.   # print('dt1={}'.format(self.dt1))
                if self.dt1 < dt_min:                    
                    break # print('dt too small')
            else:
                # if higher, continue steps along the primary
                RV_old = RV_new


#%%-----------------------------------------------------------------------------
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
        # ********************************************************* #               
        tr.pass_fan(r_aim, E, B, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z, no_intersect=True, 
            no_out_of_bounds=True, invisible_wall_x=geom.r_dict[target][0]+0.2)
        # pass trajectory to the target
        # ********************************************************* #               
        tr.pass_to_target(r_aim, E, B, geom, stop_plane_n=stop_plane_n, eps_xy=eps_xy, eps_z=eps_z, dt_min=dt_min, no_intersect=True, 
            no_out_of_bounds=True, invisible_wall_x=geom.r_dict[target][0]+0.2)
        # ********************************************************* #               

        print('IsAimXY = ', tr.IsAimXY); print('IsAimZ = ', tr.IsAimZ)
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
                print('UB2 OLD = {:.2f}, z_aim - z = {:.4f} m'.format(UB2, dz))
                
                UB2_old = UB2 
                UB2 = UB2 - dUB2*dz
                if np.isnan(UB2): 
                    tr.print_log("dUB2 = %f" % dUB2); tr.print_log("dz = %f" % dz); tr.print_log("UB2_old = %f" % UB2_old)
                    
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
    
'''    
pass_fan         calls    pass_prim    
pass_to_target   calls    pass_sec
optimize_B2      calls    pass_fan   pass_to_target
'''
#%%

def delta(r_target, r, coords): 
    if coords == 'xy':
        return np.linalg.norm(r_target[:2]-r[:2]) * np.sign(np.cross(r[:2], r_target[:2]))
    elif coords == 'z': 
        r_target[2] - r[2]
    else: 
        raise Exception("delta: invalid coord")

def optimize(tr, optimizator): 
    while True: 
        tr.pass_sec(RV0, r_target, E, B, geom, tmax=tmax, eps_xy=eps_xy, eps_z=eps_z)
        if not optimizator.change_voltages(tr): 
            return optimizator.voltage_list 

class OptimizerA3: 
    def __init__(self, r_target): 
        self.plates = 'A3'
        self.voltages_list = []
        self.r_target = r_target
        self.count = 0
        
    def change_voltages(self, tr): 
        last_r = tr.RV_sec[-1, 0:3]
        drXY = delta(r_target, last_r, 'xy') 

        if  abs(drXY) < 0.01: 
            self.voltage_list.append(tr.U['A3'])
            return True
        if count > 10000: 
            return True

        tr.U['A3'] = tr.U['A3']  + dUA3*drXY
        count += 1
        return False

        

def optimize_A3B3(tr, geom, UA3, UB3, dUA3, dUB3,
                  E, B, dt, target='slit', UA3_max=50., UB3_max=50.,
                  eps_xy=1e-3, eps_z=1e-3):
    '''
    get voltages on A3 and B3 plates to get into target
    '''
    
    #set target
    r_target = geom.r_dict[target]
    stop_plane_n = geom.plates_dict['an'].slit_plane_n
    
    #reset variables
    tr.dt1, tr.dt2 = dt, dt
    tmax = 9e-5
    tr.IsAimXY, tr.IsAimZ = False, False  # ; vltg_fail = False
    n_stepsA3 = 0
    
    RV0 = np.array([tr.RV_sec[0]])  #starting point

    while not (tr.IsAimXY and tr.IsAimZ):
        #set voltages
        tr.U['A3'], tr.U['B3'] = UA3, UB3
        #passing traj
        tr.pass_sec(RV0, r_target, E, B, geom, tmax=tmax, eps_xy=eps_xy, eps_z=eps_z)

        #calculate discrepancy from target
        drXY = delta(r_target, tr.RV_sec[-1], 'xy') 
        #set new voltage
        UA3 = UA3 + dUA3*drXY; n_stepsA3 += 1
        
        #check A3 if nan, too big or too many steps

        if abs(drXY) < 0.01:
            if tr.IntersectGeometrySec['A3']:
                return tr, True # vltg_fail = True

            n_stepsZ = 0
            while not tr.IsAimZ:

                tr.U['A3'], tr.U['B3'] = UA3, UB3
                tr.pass_sec(RV0, r_target, E, B, geom, tmax=tmax, eps_xy=eps_xy, eps_z=eps_z)
                dz = r_target[2] - tr.RV_sec[-1, 2]
                UB3 = UB3 - dUB3*dz; n_stepsZ += 1
                if ( abs(UB3) > UB3_max ) or (n_stepsZ > 100):
                    return tr, True # vltg_fail = True

            n_stepsA3 = 0
            dz =  delta(r_target, tr.RV_sec[-1], 'z') 


    return tr, False # vltg_fail = False

