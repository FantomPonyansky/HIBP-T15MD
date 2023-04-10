# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:24:00 2023

@author: Krohalev_OD
"""

# %% Optimize Secondary Beamline
t1 = time.time()
# define list of trajectories that hit slit
traj_list_a3b3 = []
traj_list_failed_to_pass_sec = []
if optimizeA3B3:
    print('\n Secondary beamline optimization')
    for tr in copy.deepcopy(traj_list_passed):
        tr, vltg_fail = hb.optimize_A3B3(tr, geomT15, UA3, UB3, dUA3, dUB3,
                                          E, B, dt, target='slit',  # 'aimA4'
                                          UA3_max=40., UB3_max=40.,
                                          eps_xy=1e-3, eps_z=1e-3)
        # check geometry intersection and voltage failure
        if not (True in tr.IntersectGeometrySec.values()) and not vltg_fail:
            traj_list_a3b3.append(tr)
            print('\n Trajectory saved')
            UA3 = tr.U['A3']
            UB3 = tr.U['B3']
        else:
            traj_list_failed_to_pass_sec.append(tr)
            print('\n NOT saved')
    t2 = time.time()
    print('\n A3 & B3 voltages optimized, t = {:.1f} s\n'.format(t2-t1))
else:
    print('\n Calculating secondary beamline')
    for tr in copy.deepcopy(traj_list_passed):
        print('\nEb = {}, UA2 = {:.2f}'.format(tr.Ebeam, tr.U['A2']))
        RV0 = np.array([tr.RV_sec[0]])
        tr.pass_sec(RV0, geomT15.r_dict['slit'], E, B, geomT15,
                    stop_plane_n=geomT15.plates_dict['an'].det_plane_n,
                    tmax=9e-5, eps_xy=eps_xy, eps_z=eps_z)
        traj_list_a3b3.append(tr)
    t2 = time.time()
    print('\n Secondary beamline calculated, t = {:.1f} s\n'.format(t2-t1))

# %%
def optimize_A3B3(tr, geom, UA3, UB3, dUA3, dUB3,
                  E, B, dt, target='slit', UA3_max=50., UB3_max=50.,
                  eps_xy=1e-3, eps_z=1e-3):
    '''
    get voltages on A3 and B3 plates to get into target
    '''
    print('\nEb = {}, UA2 = {}'.format(tr.Ebeam, tr.U['A2']))
    print('Target: ' + target)
    rs = geom.r_dict[target]
    stop_plane_n = geom.plates_dict['an'].slit_plane_n

    tr.dt1 = dt
    tr.dt2 = dt
    tmax = 9e-5
    tr.IsAimXY = False
    tr.IsAimZ = False
    RV0 = np.array([tr.RV_sec[0]])

    vltg_fail = False  # flag to track voltage failure
    n_stepsA3 = 0
    while not (tr.IsAimXY and tr.IsAimZ):
        tr.U['A3'], tr.U['B3'] = UA3, UB3
        tr.pass_sec(RV0, rs, E, B, geom,  #stop_plane_n=stop_plane_n,
                    tmax=tmax, eps_xy=eps_xy, eps_z=eps_z)

        drXY = np.linalg.norm(rs[:2]-tr.RV_sec[-1, :2]) * \
            np.sign(np.cross(tr.RV_sec[-1, :2], rs[:2]))
            
        print('\n UA3 OLD = {:.2f} kV, dr XY = {:.4f} m'.format(UA3, drXY))
        print('IsAimXY = ', tr.IsAimXY)
        # if drXY < 1e-2:
        #     dUA3 = 10.0

        UA3 = UA3 + dUA3*drXY
        print('UA3 NEW = {:.2f} kV'.format(UA3))
        n_stepsA3 += 1

        if np.isnan(UA3):
            print('ALPHA3 failed, nan value in UA3')
            vltg_fail = True
            return tr, vltg_fail
        if abs(UA3) > UA3_max:
            print('ALPHA3 failed, voltage too high')
            vltg_fail = True
            return tr, vltg_fail
        if n_stepsA3 > 100:
            print('ALPHA3 failed, too many steps')
            vltg_fail = True
            return tr, vltg_fail

        # dz = rs[2] - tr.RV_sec[-1, 2]
        # print('\n UB3 OLD = {:.2f} kV, dZ = {:.4f} m'.format(UB3, dz))
        if abs(drXY) < 0.01:
            if tr.IntersectGeometrySec['A3']:
                print('BAD A3!')
                vltg_fail = True
                return tr, vltg_fail
            n_stepsZ = 0
            while not tr.IsAimZ:
                print('pushing Z direction')
                tr.U['A3'], tr.U['B3'] = UA3, UB3
                tr.pass_sec(RV0, rs, E, B, geom, #stop_plane_n=stop_plane_n,     #!!! originally no stop plane was given
                            tmax=tmax, eps_xy=eps_xy, eps_z=eps_z)
                # tr.IsAimZ = True  # if you want to skip UB3 calculation
                dz = rs[2] - tr.RV_sec[-1, 2]
                print(' UB3 OLD = {:.2f} kV, dZ = {:.4f} m'
                      .format(UB3, dz))
                print('IsAimXY = ', tr.IsAimXY)
                print('IsAimZ = ', tr.IsAimZ)

                UB3 = UB3 - dUB3*dz
                n_stepsZ += 1
                if abs(UB3) > UB3_max:
                    print('BETA3 failed, voltage too high')
                    vltg_fail = True
                    return tr, vltg_fail
                if n_stepsZ > 100:
                    print('BETA3 failed, too many steps')
                    vltg_fail = True
                    return tr, vltg_fail
                # print('UB3 NEW = {:.2f} kV'.format(UB3))
            n_stepsA3 = 0
            print('n_stepsZ = ', n_stepsZ)
            dz = rs[2] - tr.RV_sec[-1, 2]
            print('UB3 NEW = {:.2f} kV, dZ = {:.4f} m'.format(UB3, dz))

    return tr, vltg_fail

# %%
def optimize_A4(tr, geom, UA4, dUA4, E, B, dt, eps_alpha=0.1):
    '''
    get voltages on A4 to get proper alpha angle at the entrance to analyzer
    '''
    print('\n A4 optimization')
    print('\nEb = {}, UA2 = {}'.format(tr.Ebeam, tr.U['A2']))

    rs = geom.r_dict['slit']
    stop_plane_n = geom.plates_dict['an'].slit_plane_n
    alpha_target = geom.angles_dict['an'][0]

    tr.dt1 = dt
    tr.dt2 = dt
    tmax = 9e-5
    # tr.IsAimXY = False
    # tr.IsAimZ = False
    RV0 = np.array([tr.RV_sec[0]])
    V_last = tr.RV_sec[-1][3:]
    alpha, beta = calc_angles(V_last)
    dalpha = alpha_target - alpha
    n_stepsA4 = 0
    while (abs(alpha - alpha_target) > eps_alpha):
        tr.U['A4'] = UA4
        tr.pass_sec(RV0, rs, E, B, geom,  # stop_plane_n=stop_plane_n,
                    tmax=tmax, eps_xy=1e-2, eps_z=1e-2)

        V_last = tr.RV_sec[-1][3:]
        alpha, beta = calc_angles(V_last)
        dalpha = alpha_target - alpha
        print('\n UA4 OLD = {:.2f} kV, dalpha = {:.4f} deg'.format(UA4, dalpha))
        
        #_, drXY, dz = delta_projection(tr.RV_sec[-1, :2], rs[:2], stop_plane_n)
        drXY = np.linalg.norm(rs[:2]-tr.RV_sec[-1, :2]) * np.sign(np.cross(tr.RV_sec[-1, :2], rs[:2]))
        dz = rs[2] - tr.RV_sec[-1, 2]
        print('dr XY = {:.4f} m, dz = {:.4f} m'.format(drXY, dz))

        UA4 = UA4 + dUA4*dalpha
        print('UA4 NEW = {:.2f} kV'.format(UA4))
        n_stepsA4 += 1

        if abs(UA4) > 50.:
            print('ALPHA4 failed, voltage too high')
            return tr
        if n_stepsA4 > 100:
            print('ALPHA4 failed, too many steps')
            return tr

    return tr

def optimize_sec(E, B, geom, traj_list, optimization_mode="center", 
                 target='slit', max_voltages=[40., 40., 40.], eps_xy=1e-3, 
                 eps_z=1e-3, tmax=9e-5, U0 = [0., 0., 0.], dU = [7, 10, 5.]):
    
    #create list for results
    traj_list_passed = []
    
    #create 3 optimizers
    opt_A3 = optimizers.Optimizer('zone', dU[0], max_voltages[0], target, geom, 'A3')
    opt_B3 = optimizers.Optimizer(optimization_mode, dU[1], max_voltages[1], target, geom, 'B3')
    opt_A4 = optimizers.Optimizer('center',            dU[2], max_voltages[2], target, geom, 'A4')
    
    for tr in traj_list:
        
        print('\nEb = {}, UA2 = {}'.format(tr.Ebeam, tr.U['A2']))
        print('Target: ' + target)
        
        #set start point and voltages
        RV0 = np.array([tr.RV_sec[0]])
        tr.U['A3'], tr.U['B3'], tr.U['A4'] = U0
        
        #reset optimizers
        A3_voltages, B3_voltages, A4_voltages = [], [], []
        opt_A3.voltages_list, opt_B3.voltages_list, opt_A4.voltages_list = [], [], []
        opt_A3.count, opt_B3.count, opt_A4.count = 0, 0, 0
        opt_A3.U_step, opt_B3.U_step, opt_A4.U_step = None, None, None
        
        #optimize plates one by one
        A3_voltages = optimize(tr, opt_A3, E, B, geom, RV0, tmax, eps_xy, eps_z)
        if A3_voltages:
            print("\nOptimizing B3")
            #pass A3 voltages to B3 optimizer
            opt_B3.previous_voltages['A3'] = A3_voltages
            B3_voltages = optimize(tr, opt_B3, E, B, geom, RV0, tmax, eps_xy, eps_z)
            print("B3 voltages: ", B3_voltages)
        if B3_voltages:
            print("\nOptimizing A4")
            #pass A3 and B3 voltages to A4 optimizer
            opt_A4.previous_voltages['A3'] = A3_voltages
            opt_A4.previous_voltages['B3'] = B3_voltages
            A4_voltages = optimize(tr, opt_A4, E, B, geom, RV0, tmax, eps_xy, eps_z)
        
        #finish optimization and write down results
        if (not A4_voltages) or (True in tr.IntersectGeometrySec.values()):
            print("\nOptimization failed, trajectory NOT saved\n")
        else: 
            print("\nTrajectory optimized and saved\n")
            traj_list_passed.append(tr)
            
    return traj_list_passed

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