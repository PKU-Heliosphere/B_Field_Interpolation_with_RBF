"""
@Filename: RBF_B.py
@Purpose: interpolate magnetic field vector based on a radial basis function (RBF) method, \
          which can keep the magnetic field locally divergence-free.
@Input: grid data and position of target point. 
@Output：magnetic field vector at target point in Spherical Cordinates. shape:[Bx, By, Bz]
@Coordinate: Cartesian coordinates
@More information:   
    1. The RBF interpolation procedure requires only 27 (=3x3x3) data points around the target point,\
        independent of the grid. \
       The role of the grid data in the input and \
        the subscripts of the nearest points is to find the 3x3x3 data points closest to the target point.
    2. Regardless of whether the data grid is in Cartesian or Spherical coordinates, \
        the position and magnetic field vectors used in the interpolation process \
         need to be transferred to Cartesian coordinates.
@References:
    Fortran code provided by Zhang et al. (2021)
    Zhang et al., 2021  doi: 10.3847/1538-4365/ac1e29
    Bras et al., 2019  doi: 10.1002/fld.4767
@Author: Chuanpeng Hou, Jiansen He.
@E-mail: jshept@pku.edu.cn
@Date of Last Change: 2022-02-20
"""

import numpy as np
# this is the main function
def magnetic_RBF_interp(x_data,y_data,z_data,
                       Bx_data,By_data,Bz_data,
                       xi, yi, zi,
                       x_index_of_nearest_point,
                       y_index_of_nearest_point,
                       z_index_of_nearest_point,
                       nv=27, m=26):
    '''
    :param x_data, y_data, z_data: grid point positions. 1D array
    :param Bx_data, By_data, Bz_data: magnetic data on mesh grid point. 3D array
    :param xi, yi, zi: target point position. float. example:[1.0, 2.0, 3.0]
    :param x_index_of_nearest_point: to find the 3x3x3 data points closest to the target point
    :param y_index_of_nearest_point: to find the 3x3x3 data points closest to the target point
    :param z_index_of_nearest_point: to find the 3x3x3 data points closest to the target point
    :param nv: nv = 27
    :param m: m = 26
    :return: magnetic vector on target position, 1D vector. example: [1.0, -1.0, 3.0]
    '''

    ''' 
    step1. Here's some preparation work.
           Obtain the magnetic field and position data and save them in an array of length 27.
    '''
    Bi = np.zeros(3)
    cout = 0
    B_arr = np.zeros((3, nv))
    x = np.zeros(nv)
    y = np.zeros(nv)
    z = np.zeros(nv)
    for i in range(x_index_of_nearest_point-1, x_index_of_nearest_point+2):
        for j in range(y_index_of_nearest_point - 1, y_index_of_nearest_point + 2):
            for k in range(z_index_of_nearest_point - 1, z_index_of_nearest_point + 2):
                # print(cout)
                B_arr[0, cout] = Bx_data[i, j, k]
                B_arr[1, cout] = By_data[i, j, k]
                B_arr[2, cout] = Bz_data[i, j, k]
                
                x[cout] = x_data[i]
                y[cout] = y_data[j]
                z[cout] = z_data[k]
                cout = cout + 1
            # end
        # end
    # end

    '''
    Step2. Starting Magnetic field interpolation
    '''
    # Do not change the following:
    r = np.zeros((nv))

    for j in range(nv):
        r[j] = np.sqrt((xi-x[j])**2 + (yi-y[j])**2 + (zi-z[j])**2)
    r_max = np.nanmax(r)
    r_support = 2 * r_max

    PHI = np.zeros((3,3,nv,nv))
    for j in range(nv):
        for k in range(nv):
            psi = RBF_function_partial1(r_support, x[j], y[j], z[j],
                                        x[k], y[k], z[k])
            for j1 in [0, 1, 2]:
                for k1 in [0, 1, 2]:
                    PHI[j1, k1, j, k] = psi[j1, k1]
    vector = np.zeros(3)
    P = np.zeros((3*nv, m))
    for j in range(nv):
        for k in range(m):
            psi = polynomial_magnetic(x[j], y[j], z[j], k)
            for j1 in [0, 1, 2]:
                vector[j1] = psi[j1]
            j1 = (j)*3
            # print(vector)
            # input()
            # P[j1+1:j1+4,k] = vector[:]
            P[j1 + 0, k] = vector[0]
            P[j1 + 1, k] = vector[1]
            P[j1 + 2, k] = vector[2]


    PT = np.transpose(P)
    A = np.zeros((3 * nv + m, 3 * nv + m))
    c = np.zeros(3 * nv + m)
    for j in range(0, 3*nv-2, 3):
        for k in range(0, 3 * nv - 2, 3):
            for j1 in [j, j+1, j+2]:
                for k1 in [k, k + 1, k + 2]:
                    # print(k1)
                    A[j1, k1] = PHI[j1-j, k1-k, int(j/3), int(k/3)] # 不需要加1

    for j in range(3*nv):
        for k in range(3*nv, 3*nv+m):
            A[j, k] = P[j, k-3*nv]
    for j in range(3*nv, 3*nv+m):
        for k in range(3*nv):
            A[j, k] = PT[j-3*nv, k]
    for j in range(0, 3*nv-2, 3):
        for j1 in [j, j+1, j+2]:
            c[j1] = B_arr[j1-j, int(j/3)] # 不需要加1
    w = matrix_solver(A, c)
    for j in range(nv):
        psi = RBF_function_partial1(r_support, xi, yi, zi,
                                    x[j], y[j], z[j])
        for k in range(3):
            Bi[0] = Bi[0] + psi[0, k] * w[(j) * 3 + k]
            Bi[1] = Bi[1] + psi[1, k] * w[(j) * 3 + k]
            Bi[2] = Bi[2] + psi[2, k] * w[(j) * 3 + k]
    for j in range(m):
        psi = polynomial_magnetic(xi, yi, zi, j)
        Bi[0] = Bi[0] + psi[0] * w[3 * nv + j]
        Bi[1] = Bi[1] + psi[1] * w[3 * nv + j]
        Bi[2] = Bi[2] + psi[2] * w[3 * nv + j]
    return Bi  # interpolation result [Bx,By,Bz]

# some functions used in the main function
def RBF_function_partial1(r_support, x,y,z,xj,yj,zj):
    """
    :param r_support: single value which is two times of maximum value of distance between target point and neighbor points.
    :param x: x of neighbor points in Cartesian Coordinates. shape: 1D
    :param y: y of neighbor points in Cartesian Coordinates. shape: 1D
    :param z: z of neighbor points in Cartesian Coordinates. shape: 1D
    :param xj: x of target point in Cartesian Coordinates. single value
    :param yj: y of target point in Cartesian Coordinates. single value
    :param zj: z of target point in Cartesian Coordinates. single value
    :return: partial: the basis function
    """
    alph = 1.0 / r_support ** 2
    temp = np.exp(-alph * ((x - xj) ** 2 + (y - yj) ** 2 + (z - zj) ** 2))
    pxx = temp * (4 * alph ** 2 * (x - xj) ** 2 - 2 * alph)
    pyy = temp * (4 * alph ** 2 * (y - yj) ** 2 - 2 * alph)
    pzz = temp * (4 * alph ** 2 * (z - zj) ** 2 - 2 * alph)
    pxy = 4 * alph ** 2 * (x - xj) * (y - yj) * temp
    pxz = 4 * alph ** 2 * (x - xj) * (z - zj) * temp
    pyz = 4 * alph ** 2 * (y - yj) * (z - zj) * temp

    partial = np.zeros([3, 3])
    partial[0, 0] = - pyy - pzz
    partial[0, 1] = pxy
    partial[0, 2] = pxz
    partial[1, 0] = pxy
    partial[1, 1] = -pxx - pzz
    partial[1, 2] = pyz
    partial[2, 0] = pxz
    partial[2, 1] = pyz
    partial[2, 2] = - pxx -pyy
    return partial
def polynomial_magnetic(x,y,z,m):
    """
    :param x:  x of neighbor points in Cartesian Coordinates. shape: 1D
    :param y:  y of neighbor points in Cartesian Coordinates. shape: 1D
    :param z:  z of neighbor points in Cartesian Coordinates. shape: 1D
    :param m:  the number of polynomial functions. sigle value
    :return: p: polynomial functions
    """
    p = np.zeros(3)
    if m == 0:
        p[0] = 1
    elif m == 1:
        p[1] = 1
    elif m == 2:
        p[2] = 1
    elif m == 3:
        p[0] = y
    elif m == 4:
        p[0] = z
    elif m == 5:
        p[1] = x
    elif m == 6:
        p[1] = z
    elif m == 7:
        p[2] = x
    elif m == 8:
        p[2] = y
    elif m == 9:
        p[0] = x
        p[1] = -y
    elif m == 10:
        p[0] = x
        p[2] = -z
    elif m == 11:
        p[2] = x**2
    elif m == 12:
        p[1] = x**2
    elif m == 13:
        p[2] = y**2
    elif m == 14:
        p[0] = y**2
    elif m == 15:
        p[1] = z**2
    elif m == 16:
        p[0] = z**2
    elif m == 17:
        p[2] = x*y
    elif m == 18:
        p[0] = x**2
        p[1] = -2*x*y
    elif m == 19:
        p[0] = -x**2
        p[2] = 2*x*z
    elif m == 20:
        p[1] = x*z
    elif m == 21:
        p[0] = 2*x*y
        p[1] = -y**2
    elif m == 22:
        p[0] = -2*x*z
        p[2] = z**2
    elif m == 23:
        p[1] = y**2
        p[2] = -2*y*z
    elif m == 24:
        p[0] = y*z
    elif m == 25:
        p[1] = 2*y*z
        p[2] = -z**2
    return p
def matrix_solver(A, c):
    """
    :param A: matrix
    :param c: vector
    :return: w: =A^(-1)c, a solution of linear system of equations
    """
    # A*w=c -> w = A^(-1) * c
    A_inv = np.linalg.inv(A)
    w = np.dot(A_inv, c)
    return w
