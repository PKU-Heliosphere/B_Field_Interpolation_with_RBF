# put RBF_B.py and Verify_RBF_B_function.py into one directory, \
# and just run Verify_RBF_B_function.py directly
'''
@Filename: Verify_RBF_B_function.py
@Purpose: Test RBF interpolation code for magnetic field distribution with magnetic dipoles
@Output：Comparison figure of theoretical values and interpolation results.
@Author: Chuanpeng Hou, Jiansen He.
@E-mail: jshept@pku.edu.cn
@Date of Last Change: 2022-02-20
'''
import numpy as np
import matplotlib.pyplot as plt
from RBF_B import magnetic_RBF_interp

def delta_function(r_norm):
    if r_norm<1e-10:
        return 1
    else:
        return 0
def get_B(r):
    mu0 = 4*np.pi*1e-7
    m = np.array([0, 0, 0.5]) *1e12
    # r = np.array([1, 3, 15])
    r_norm = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    m_dot_r = np.dot(m, r)/r_norm
    B = mu0/4/np.pi / r_norm**3 * ((3 * m_dot_r) * r / r_norm - m) + 2*mu0 * m /3 * delta_function(r_norm)**3
    return B
def get_xyz(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    xyz = np.array([x, y, z])
    return xyz

r_arr = np.linspace(40,100,5)
theta_arr = np.linspace(0, np.pi, 10)
phi_arr = np.linspace(0, 2*np.pi, 1)
x_grid = np.linspace(10, 90, 10)
y_grid = np.linspace(40, 90, 9)
z_grid = np.linspace(40, 90, 10)
# print(x_grid)
# y_grid = y_grid *0
xv, yv, zv = np.meshgrid(x_grid, y_grid, z_grid)
# phi_arr = np.array([0,np.pi])
# B = np.zeros((len(r_arr)*len(theta_arr)*len(phi_arr), 3))
# Bx_arr = np.zeros((len(r_arr)*len(theta_arr)*len(phi_arr)))
# By_arr = np.zeros((len(r_arr)*len(theta_arr)*len(phi_arr)))
# Bz_arr = np.zeros((len(r_arr)*len(theta_arr)*len(phi_arr)))

x_arr = np.zeros((len(x_grid)*len(y_grid)*len(z_grid)))
y_arr = np.zeros((len(x_grid)*len(y_grid)*len(z_grid)))
z_arr = np.zeros((len(x_grid)*len(y_grid)*len(z_grid)))

B = np.zeros((len(x_grid)*len(y_grid)*len(z_grid), 3))
Bx_arr = np.zeros((len(x_grid)*len(y_grid)*len(z_grid)))
By_arr = np.zeros((len(x_grid)*len(y_grid)*len(z_grid)))
Bz_arr = np.zeros((len(x_grid)*len(y_grid)*len(z_grid)))

Bx_grid_T = np.zeros((len(x_grid),len(y_grid),len(z_grid)))
By_grid_T = np.zeros((len(x_grid),len(y_grid),len(z_grid)))
Bz_grid_T = np.zeros((len(x_grid),len(y_grid),len(z_grid)))
Bx_grid = xv*0
By_grid = xv*0
Bz_grid = xv*0
# print()
# x_arr = np.zeros((len(r_arr)*len(theta_arr)*len(phi_arr)))
# y_arr = np.zeros((len(r_arr)*len(theta_arr)*len(phi_arr)))
# z_arr = np.zeros((len(r_arr)*len(theta_arr)*len(phi_arr)))
cout = 0
# ax = plt.figure().add_subplot(projection='3d')
# for r in r_arr:
#     for theta in theta_arr:
#         for phi in phi_arr:

shape_xv = np.shape(xv)
xx_len = shape_xv[0]
yy_len = shape_xv[1]
zz_len = shape_xv[2]

for ii in range(xx_len):
    for jj in range(yy_len):
        for kk in range(zz_len):
            # xyz = get_xyz(r, theta, phi)
            xyz = np.array([xv[ii,jj,kk], yv[ii,jj,kk], zv[ii,jj,kk]])
            # print(xyz)
            B = get_B(xyz)
            x_arr[cout] = xyz[0]
            y_arr[cout] = xyz[1]
            z_arr[cout] = xyz[2]
            Bx_arr[cout] = B[0]
            By_arr[cout] = B[1]
            Bz_arr[cout] = B[2]
            Bx_grid[ii, jj, kk] = B[0]
            By_grid[ii, jj, kk] = B[1]
            Bz_grid[ii, jj, kk] = B[2]
            cout = cout + 1
            # ax.quiver(xyz[0], xyz[1], xyz[2], B[cout-1,0], B[cout-1,1], B[cout-1,2],length=0.5, normalize=True)
            plt.quiver(xyz[0], xyz[2], Bx_arr[cout - 1], Bz_arr[cout - 1])
plt.xlabel('x')
plt.ylabel('y')
# plt.show()
num = 0
for x in np.linspace(25,70,15):
    ri =np.array([x,55,54])
    x_index_of_nearest_point = np.argmin(np.abs(x_grid-ri[0]))
    y_index_of_nearest_point = np.argmin(np.abs(y_grid-ri[1]))
    z_index_of_nearest_point = np.argmin(np.abs(z_grid-ri[2]))
    B_true = get_B(ri)
    nv = 27
    m = 26
    test_pos = np.array([x_grid[x_index_of_nearest_point],y_grid[y_index_of_nearest_point],z_grid[z_index_of_nearest_point]])
    test_Bx = Bx_grid[y_index_of_nearest_point,x_index_of_nearest_point,z_index_of_nearest_point]
    test_By = By_grid[y_index_of_nearest_point,x_index_of_nearest_point,z_index_of_nearest_point]
    test_Bz = Bz_grid[y_index_of_nearest_point,x_index_of_nearest_point,z_index_of_nearest_point]
    B= get_B(test_pos)

    # print(test_pos,B)
    # print('test_B: ',[test_Bx,test_By,test_Bz])

    for i in range(zz_len):
        Bx_grid_T[:,:,i] = np.transpose(Bx_grid[:,:,i])
        By_grid_T[:,:,i] = np.transpose(By_grid[:,:,i])
        Bz_grid_T[:,:,i] = np.transpose(Bz_grid[:,:,i])
    # input()
    B = magnetic_RBF_interp(x_grid,y_grid,z_grid,
                       Bx_grid_T,By_grid_T,Bz_grid_T,
                       ri[0], ri[1], ri[2],
                           x_index_of_nearest_point,
                           y_index_of_nearest_point,
                           z_index_of_nearest_point,
                           nv, m)
    plt.subplot(3,1,1)

    plt.scatter(num,B[0],c='red', alpha=1,marker='*')
    plt.scatter(num, B_true[0], c='black', marker='o', alpha=0.4)
    plt.xlabel('sampling #')
    plt.ylabel('Bx')
    plt.legend(("interpolation", "theory"), loc="upper right")

    plt.subplot(3,1,2)
    plt.scatter(num,B[1],c='red',marker='*')
    plt.scatter(num, B_true[1], c='black', marker='o', alpha=0.4)
    plt.xlabel('sampling #')
    plt.ylabel('By')
    plt.legend(("interpolation", "theory"), loc="upper right")

    plt.subplot(3,1,3)
    plt.scatter(num,B[2],c='red', alpha=1,marker='*')
    plt.scatter(num, B_true[2], c='black', marker='o', alpha=0.4)
    plt.xlabel('sampling #')
    plt.ylabel('Bz')
    num = num +1
    plt.legend(("interpolation", "theory"), loc="upper right")
plt.show()
print('B_theory: ', B_true)
print('B_interp: ', B)




