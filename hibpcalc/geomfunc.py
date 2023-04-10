# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 15:33:36 2023

@author: Krohalev_OD
"""
#%% imports
import numpy as np
import matplotlib.pyplot as plt

import numba 

# %%
def calc_normal(point1, point2, point3):
    '''
    calculate vector normal to a plane defined by 3 points
    '''
    plane_n = np.cross(point1 - point2, point1 - point3)
    return plane_n/np.linalg.norm(plane_n)

#%% auxilliary functions
def front_rect(edges): 
    return [ edges[0, 0], edges[1, 0], edges[1, 1], edges[0, 1], ]

def rect_center(rect): 
    return 0.25*(rect[0] + rect[1] + rect[2] + rect[3] )
    
def inner_rect(rect, k=0.5): 
    k = 1.0 - k # 1: original rect, 0: center point
    center = rect_center(rect)
    return [ rect[i] + k*(center - rect[i]) for i in range(4) ]

def shifted_point_of_rect(rect, kx, ky): 
    pass
    
def plot_rect(rect, fig=None, tangage=False, **kwargs): 
    if fig is not None:
        plt.figure(fig)
    xx = [rect[i][0] for i in range(4)] 
    yy = [rect[i][1] for i in range(4)] 
    xx.append(xx[0])
    yy.append(yy[0])
    plt.plot(xx, yy, **kwargs)
    
    if tangage: 
        #p1 = 0.5*(rect[0] + rect[1])
        #p2 = 0.5*(rect[2] + rect[3])
        p1 = 0.5*(rect[1] + rect[2])
        p2 = 0.5*(rect[0] + rect[3])
        
        xx = [p1[0], p2[0]] 
        yy = [p1[1], p2[1]] 
        plt.plot(xx, yy, **kwargs)

def plot_point(pt, fig=None, ms=4, **kwargs): 
    if fig is not None:
        plt.figure(fig)
    xx = [pt[0]] 
    yy = [pt[1]] 
    plt.plot(xx, yy, 'o', ms=ms, **kwargs) # markersize
    
    

def plot_segm(p0, p1, fig=None):     
    if fig is not None:
        plt.figure(fig)
    plt.plot([  p0[0], p1[0]  ], 
             [  p0[1], p1[1]  ] )

# %%
@numba.njit()
def calc_vector(length, alpha, beta):
    '''
    calculate vector based on its length and angles
    alpha is the angle with XZ plane
    beta is the angle of rotation around Y axis
    '''
    drad = np.pi/180.  # converts degrees to radians
    x = np.cos(alpha*drad) * np.cos(beta*drad)
    y = np.sin(alpha*drad)
    z = -np.cos(alpha*drad) * np.sin(beta*drad)
    return np.array([x, y, z]) * length

@numba.njit()
def calc_angles(vector):
    '''
    calculate alpha and beta angles based on vector coords
    '''
    drad = np.pi/180.  # converts degrees to radians
    x, y, z = vector / np.linalg.norm(vector)
    # alpha = np.arcsin(y)  # rad
    alpha = np.arccos(x)  # rad
    if abs(y) > 1e-9:
        beta = np.arcsin(-np.tan(alpha) * z / y)  # rad
    elif abs(z) < 1e-9:
        beta = 0.
    elif abs(x) > 1e-9:
        beta = np.arctan(-z / x)  # rad
    else:
        beta = -np.sign(z) * np.pi/2
    return alpha/drad, beta/drad  # degrees

# %% get axes index
def get_index(axes):
    axes_dict = {'XY': (0, 1), 'XZ': (0, 2), 'ZY': (2, 1)}
    return axes_dict[axes]

# %%
@numba.njit()
def translate(input_array, xyz):
    '''
    move the vector in space
    xyz : 3 component vector
    '''
    if input_array is not None:
        input_array += np.array(xyz)

    return input_array

@numba.njit()
def rot_mx(axis=(1, 0, 0), deg=0):
    '''
    function calculates rotation matrix
    '''
    n = axis
    ca = np.cos(np.radians(deg))
    sa = np.sin(np.radians(deg))
    R = np.array([[n[0]**2*(1-ca)+ca, n[0]*n[1]*(1-ca)-n[2]*sa,
                   n[0]*n[2]*(1-ca)+n[1]*sa],

                  [n[1]*n[0]*(1-ca)+n[2]*sa, n[1]**2*(1-ca)+ca,
                   n[1]*n[2]*(1-ca)-n[0]*sa],

                  [n[2]*n[0]*(1-ca)-n[1]*sa, n[2]*n[1]*(1-ca)+n[0]*sa,
                   n[2]**2*(1-ca)+ca]])
    return R


@numba.njit()
def rotate(input_array, axis=(1, 0, 0), deg=0.):
    '''
    rotate vector around given axis by deg [degrees]
    axis : axis of rotation
    deg : angle in degrees
    '''
    if input_array is not None:
        R = rot_mx(axis, deg)
        input_array = np.dot(input_array, R.T)
    return input_array


@numba.njit()
def rotate3(input_array, plates_angles, beamline_angles, inverse=False):
    '''
    rotate vector in 3 dimentions
    plates_angles : angles of the plates
    beamline_angles : angles of the beamline axis, rotation on gamma angle
    '''
    alpha, beta, gamma = plates_angles
    beamline_axis = calc_vector(1, beamline_angles[0], beamline_angles[1])

    if inverse:
        rotated_array = rotate(input_array, axis=beamline_axis, deg=-gamma)
        rotated_array = rotate(rotated_array, axis=(0, 1, 0), deg=-beta)
        rotated_array = rotate(rotated_array, axis=(0, 0, 1), deg=-alpha)
    else:
        rotated_array = rotate(input_array, axis=(0, 0, 1), deg=alpha)
        rotated_array = rotate(rotated_array, axis=(0, 1, 0), deg=beta)
        rotated_array = rotate(rotated_array, axis=beamline_axis, deg=gamma)
    return rotated_array

def get_rotate3_mx(plates_angles, beamline_angles, inverse=False): 
    alpha, beta, gamma = plates_angles
    beamline_axis = calc_vector(1.0, beamline_angles[0], beamline_angles[1])
    
    mx1 = rotateMx( vec3D(0.0, 0.0, 1.0), np.radians(alpha) )
    mx2 = rotateMx( vec3D(0.0, 1.0, 0.0), np.radians(beta)  )
    mx3 = rotateMx( beamline_axis,        np.radians(gamma) )
    
    mx = mx3.dot( mx2.dot(mx1) )  # mx1*mx2*mx3
    return np.matrix( mx ).I.A if inverse else mx 


# %% Intersection check functions
@numba.njit()
def line_plane_intersect(planeNormal, planePoint, rayDirection,
                         rayPoint, eps=1e-6):
    '''
    function returns intersection point between plane and ray
    '''
    ndotu = np.dot(planeNormal, rayDirection)
    if abs(ndotu) < eps:
        # print('no intersection or line is within plane')
        return np.full_like(planeNormal, np.nan)
    else:
        w = rayPoint - planePoint
        si = -np.dot(planeNormal, w) / ndotu
        Psi = w + si * rayDirection + planePoint
        return Psi


@numba.njit()
def is_between(A, B, C, eps=1e-4):
    '''
    function returns True if point C is on the segment AB (between A and B)
    '''
    if np.isnan(C).any():
        return False
    # check if the points are on the same line
    crossprod = np.cross(B-A, C-A)
    if np.linalg.norm(crossprod) > eps:
        return False
    # check if the point is between
    dotprod = np.dot(B-A, C-A)
    if dotprod < 0 or dotprod > np.linalg.norm(B-A)**2:
        return False
    return True

#%% from reometry

def vec3D(x, y, z):
    return np.array([x, y, z], dtype=np.float64)

def identMx():
    return np.array([[ 1.0,  0.0,  0.0],
                      [ 0.0,  1.0,  0.0],
                      [ 0.0,  0.0,  1.0]],   dtype=np.float64)
    
def rotateMx(axis, ang):    
    res = identMx()
    s = np.sin(ang)
    c = np.cos(ang)
    x, y, z = -axis/np.linalg.norm(axis)   
    xx, yy, zz = x*x, y*y, z*z
    xy, xz, yz = x*y, x*z, y*z

    res[0, 0] = xx + (1 - xx)*c
    res[0, 1] = xy*(1 - c) + z*s
    res[0, 2] = xz*(1 - c) - y*s

    res[1, 0] = xy*(1 - c) - z*s
    res[1, 1] = yy + (1 - yy)*c
    res[1, 2] = yz*(1 - c) + x*s

    res[2, 0] = xz*(1 - c) + y*s
    res[2, 1] = yz*(1 - c) - x*s
    res[2, 2] = zz + (1 - zz)*c

    return res

def xRotateMx(ang):
    s = np.sin(ang)
    c = np.cos(ang)
    return np.array([[ 1.0,  0.0,  0.0],
                      [ 0.0,    c,   -s],
                      [ 0.0,    s,    c]],   dtype=np.float64)

def yRotateMx(ang):
    s = np.sin(ang)
    c = np.cos(ang)    
    return np.array([[ c,    0.0,    s],
                      [ 0.0,  1.0,  0.0],
                      [-s,    0.0,    c]],   dtype=np.float64)

def zRotateMx(ang):
    s = np.sin(ang)
    c = np.cos(ang)
    return np.array([[ c,   -s,    0.0],
                      [ s,    c,    0.0],
                      [ 0.0,  0.0,  1.0]],   dtype=np.float64)

def transformPt(pt, mx):   # 2D, 3D
    return mx.dot(pt).A1

def zwardRotateMx(vec):    # vec -> (0, 0, L)
# to transform any plane to XY plane
    v = vec
    a1 = np.arctan2(v[0], v[2])
    mx1 = yRotateMx(-a1)

    v = transformPt(v, mx1)
    a2 = np.arctan2(v[1], v[2])
    mx2 = xRotateMx(a2)  

    return mx2.dot(mx1)


def maxY_IdxAndVal(P, i, j, Y):
    if P[i][Y] > P[j][Y] :   return (i, P[i][Y])
    else:                    return (j, P[j][Y])

def minY_onlyVal(P, i, j, Y):
    if P[i][Y] < P[j][Y] :   return P[i][Y]
    else:                    return P[j][Y]

def ptInPolygon2D(P, pt, X=0, Y=1): 
    '''
    P - point list

    '''
    count = len(P)
    intersect = 0
    #X = 0
    #Y = 1    
    for i in range(count): 
      j = (i+1) % count
      if ( not(    ( pt[Y] < P[i][Y] )and( pt[Y] < P[j][Y] )   )and
           not(    ( pt[Y] > P[i][Y] )and( pt[Y] > P[j][Y] )   )and
           not(      P[i][Y] ==  P[j][Y]                       )    ):
        (k, maxY) = maxY_IdxAndVal(P, i, j)
        
        if maxY == pt[Y]: 
            if P[k][X] > pt[X]: 
               intersect += 1 
        else:  
          # if (not ( min(P[i][Y], P[j][Y]) == pt[Y] )): # ??? why not? 
            if (not ( minY_onlyVal(P, i, j) == pt[Y] )): 
               t = (pt[Y] - P[i][Y] )/( P[j][Y] - P[i][Y] )
               if (  (t > 0.0)and
                     (t < 1.0)and
                     (P[i][X] + t*(P[j][X]-P[i][X]) > pt[X])  ): 
                   intersect +=1
    
    return intersect % 2 == 1 


def _ptInPolygon3D_(P, pt): #  ???
    '''no check if point is in the plane of polygon 
    This function is purposed only for the points that 
    lie on the plane of the polygon
    '''
    normal = np.cross( P[1] - P[0], P[-1] - P[0])
    k = np.argmax(normal)
    X, Y = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[k]
    return ptInPolygon2D(P, pt, X, Y)

#%%

@numba.njit()
def plane_segment_intersect(planeNormal, planePoint, segmPoint0, segmPoint1, eps=1e-6): #!!! to be tested
    '''
    function returns intersection segment between plane and ray
    '''
    segmVector = segmPoint1 - segmPoint0
    planeC = - planeNormal.dot(planePoint)
    
    up = planeNormal.dot(segmPoint0) + planeC
    dn = planeNormal.dot(segmVector)
    if abs(dn) < eps: 
        # return np.full_like(segmPoint0, np.nan) 
        return None

    t = - up/dn
    intersectLinePlane = segmVector*t + segmPoint0
    
    if 0.0 <= t <= 1.0: 
        return intersectLinePlane
    else: 
        return None

def _plane_ray_intersect(planeNormal, planePoint, rayPoint, rayVector, eps=1e-6): #!!! to be tested
    '''
    function returns intersection segment between plane and ray
    '''
    planeC = - planeNormal.dot(planePoint)    
    up = planeNormal.dot(rayPoint) + planeC
    dn = planeNormal.dot(rayVector)
    if abs(dn) < eps:  # parallel 
        return None, None

    t = - up/dn
    intersectLinePlane = rayVector*t + rayPoint
    return intersectLinePlane, t

@numba.njit()
def order(A, B, C):
    '''
    if counterclockwise return True
    if clockwise return False
    '''
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])


@numba.njit()
def is_intersect(A, B, C, D):  # doesn't work with collinear case
    '''
    function returns true if line segments AB and CD intersect
    '''
    # Return true if line segments AB and CD intersect
    return order(A, C, D) != order(B, C, D) and \
        order(A, B, C) != order(A, B, D)


@numba.njit()
def segm_intersect(A, B, C, D):
    '''
    function calculates intersection point between vectors AB and CD
    in case AB and CD intersect
    '''
    # define vectors
    AB, CA, CD = B - A, A - C, D - C
    return A + AB * (np.cross(CD, CA) / np.cross(AB, CD))


@numba.njit()
def segm_poly_intersect(polygon_coords, segment_coords):
    '''
    check segment and polygon intersection
    '''
    polygon_normal = np.cross(polygon_coords[2, 0:3]-polygon_coords[0, 0:3],
                              polygon_coords[1, 0:3]-polygon_coords[0, 0:3])
    polygon_normal = polygon_normal/np.linalg.norm(polygon_normal)
    # find the intersection point between polygon plane and segment line
    intersect_coords = line_plane_intersect(polygon_normal,
                                            polygon_coords[2, 0:3],
                                            segment_coords[1, 0:3] -
                                            segment_coords[0, 0:3],
                                            segment_coords[0, 0:3])
    if np.isnan(intersect_coords).any():
        return False
    if not is_between(segment_coords[0, 0:3], segment_coords[1, 0:3],
                      intersect_coords):
        return False
    # go to 2D, exclude the maximum coordinate
    i = np.argmax(np.abs(polygon_normal))
    inds = np.array([0, 1, 2])  # indexes for 3d corrds
    inds_flat = np.where(inds != i)[0]
    polygon_coords_flat = polygon_coords[:, inds_flat]
    intersect_coords_flat = intersect_coords[inds_flat]
    # define a rectange which contains the flat poly
    xmin = np.min(polygon_coords_flat[:, 0])
    xmax = np.max(polygon_coords_flat[:, 0])
    ymin = np.min(polygon_coords_flat[:, 1])
    ymax = np.max(polygon_coords_flat[:, 1])
    xi, yi = intersect_coords_flat
    # simple check if a point is inside a rectangle
    if (xi < xmin or xi > xmax or yi < ymin or yi > ymax):
        return False
    # ray casting algorithm
    # set up a point outside the flat poly
    point_out = np.array([xmin - 0.01, ymin - 0.01])
    # calculate the number of intersections between ray and the poly sides
    intersections = 0
    for i in range(polygon_coords_flat.shape[0]):
        if is_intersect(point_out, intersect_coords_flat,
                        polygon_coords_flat[i-1], polygon_coords_flat[i]):
            intersections += 1
    # if the number of intersections is odd then the point is inside
    if intersections % 2 == 0:
        return False
    else:
        return True

    # p = path.Path(polygon_coords_flat)
    # return p.contains_point(intersect_coords_flat)

    # check projections on XY and XZ planes
    # pXY = path.Path(polygon_coords[:, [0, 1]])  # XY plane
    # pXZ = path.Path(polygon_coords[:, [0, 2]])  # XZ plane
    # return pXY.contains_point(intersect_coords[[0, 1]]) and \
    #     pXZ.contains_point(intersect_coords[[0, 2]]) and \
    #         is_between(segment_coords[0, 0:3], segment_coords[1, 0:3],
    #                   intersect_coords)

#%% 
def sign_eps(x, eps):
    if abs(x) < eps: 
        return 0.0
    else: 
        return np.sign(x)
    
#%% additional functions

def vert_horz_basis(plane_n): 
    eq_plane_n = np.array([0.0, 1.0, 0.0])
    vert_plane_n = np.array([0.0, 0.0, 1.0])
    
    vert_vector = np.cross( vert_plane_n, plane_n )
    horz_vector = np.cross( eq_plane_n, plane_n)
    
    return (vert_vector/np.linalg.norm(vert_vector), 
            horz_vector/np.linalg.norm(horz_vector) )

#%%


def _rect(r, plane_n, size): 
    rect = [np.array([-size, -size, 0.0]), 
            np.array([-size,  size, 0.0]), 
            np.array([ size,  size, 0.0]), 
            np.array([ size, -size, 0.0])    ]
    
    mx = zwardRotateMx(plane_n).I
    rect = [mx.dot(_r).A1 for _r in rect]
    rect = [_r + r for _r in rect]
    return rect


def plot_plane(plane_n, pt, size): 
    r = _rect(pt, plane_n, size)
    plot_rect(r)
    r = _rect(pt, plane_n, size*2)
    plot_rect(r)
    r = _rect(pt, plane_n, size*3)
    plot_rect(r)
    

