import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import RegularGridInterpolator

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

#########################################################################################

Earth_Radius = 6378

#########################################################################################

def data_for_cylinder_along_z(center_x,center_y,radius,height_z, earth_radius):
    z = np.linspace(earth_radius-height_z, earth_radius, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

#########################################################################################

def Density(R,h,rho):
    dens = np.zeros(R.shape);
    x = R/Earth_Radius
    dens = np.where(R<=1221.5, 13.0885 - 8.8381*x**2, dens)
    dens = np.where((1221.5<R) & (R<=3480), 12.5815 - x*(1.2638 + x*(3.6426 + x*5.5281)), dens)
    dens = np.where((3480<R) & (R<=5701), 7.9565 - x*(6.4761 - x*(5.5283 - x*3.0807)), dens)
    dens = np.where((5701<R) & (R<=5771), 5.3197 - 1.4836*x, dens)
    dens = np.where((5771<R) & (R<=5971), 11.2494 - 8.0298*x, dens)
    dens = np.where((5971<R) & (R<=6151), 7.1089 - 3.8045*x, dens)
    dens = np.where((6151<R) & (R<=6346.6), 2.691 + 0.6924*x, dens)
    dens = np.where((6346.6<R) & (R<=6356), 2.9, dens)
    dens = np.where((6356<R) & (R<=Earth_Radius-h), 2.6, dens)
    dens = np.where((Earth_Radius-h<R) & (R<=Earth_Radius), rho, dens)
    return dens

#########################################################################################

theta_EEs = np.logspace(np.log10(0.01), np.log10(90),1000)
L_tot = np.array([Earth_Radius*np.sin(2*theta_EE)/np.sin(np.pi/2-theta_EE) for theta_EE in theta_EEs])
L_ratio = np.linspace(0, 1, 10000)
Ls = L_tot[:, None]*L_ratio[None, :]

rs = np.sqrt(Ls**2+Earth_Radius**2-2*Earth_Radius*Ls*np.cos(np.pi/2-theta_EEs[:, None]))
density = Density(rs, 4, 0.92)
Xs = integrate.cumtrapz(density, Ls*1e5, axis = -1, initial = 0)

depth_interpolator = RegularGridInterpolator((theta_EEs, L_ratio), Xs)

#########################################################################################

energy = 21

Earth_Radius = 6378

detector_radius = 15
detector_depth = 2.8

angular_cut = 85

desired_events = 100000
thrown_events = 0
number_of_samples = 4000000

#events = np.empty((6,0), float)
events = np.empty((12,0), float)

#########################################################################################
#Randomly sam[ple two points on Earth's surface isotropically and draw a path between them

while events.shape[1]<desired_events:
    print(events.shape)
    thrown_events+=number_of_samples
    
    phi_local_1 = np.arccos(1-2*np.random.uniform(low = 0, high = 1, size = number_of_samples))
    theta_local_1 = np.random.uniform(size = number_of_samples)*2*np.pi
    
    phi_local_2 = np.arccos(1-2*np.random.uniform(low = 0, high = 1, size = number_of_samples))
    theta_local_2 = np.random.uniform(size = number_of_samples)*2*np.pi
    
    Earth_pos_1 = Earth_Radius*np.array([np.cos(theta_local_1)*np.sin(phi_local_1), np.sin(theta_local_1)*np.sin(phi_local_1), np.cos(phi_local_1)])
    Earth_pos_2 = Earth_Radius*np.array([np.cos(theta_local_2)*np.sin(phi_local_2), np.sin(theta_local_2)*np.sin(phi_local_2), np.cos(phi_local_2)])

#########################################################################################
#Project trajectory on x-y plane and filter only those trajectories which pass through cylinder in x-y plane
    
    slope = (Earth_pos_2[1]-Earth_pos_1[1])/(Earth_pos_2[0]-Earth_pos_1[0])
    intercept = Earth_pos_1[1]-slope*Earth_pos_1[0]
    
    a = (1+slope**2)
    b = 2*slope*intercept
    c = intercept**2-detector_radius**2
    
    root_factor = b**2-4*a*c
    within_circle = root_factor>0
    
    Earth_pos_1 = Earth_pos_1[:,within_circle]
    Earth_pos_2 = Earth_pos_2[:,within_circle]

#########################################################################################
#Calculate the points where the trajectory intersects an infinite cylinder

    a = a[within_circle]
    b = b[within_circle]
    c = c[within_circle]
    root_factor = root_factor[within_circle]
    
    slope = slope[within_circle]
    intercept = intercept[within_circle]
    
    x1, x2 = (-b-np.sqrt(root_factor))/(2*a), (-b+np.sqrt(root_factor))/(2*a)
    y1, y2 = slope*x1+intercept, slope*x2+intercept
    
    slope = (Earth_pos_2[2]-Earth_pos_1[2])/(Earth_pos_2[0]-Earth_pos_1[0])
    intercept = Earth_pos_2[2]-slope*Earth_pos_2[0]
    
    z1, z2 = slope*x1+intercept, slope*x2+intercept
    
#########################################################################################
#Filter trajectories which intersect the smaller cylinder
    
    within_cylinder_1 = np.logical_and(z1<=Earth_Radius, z2>=Earth_Radius-detector_depth)
    within_cylinder_2 = np.logical_and(z2<=Earth_Radius, z1>=Earth_Radius-detector_depth)
    
    within_cylinder = np.logical_or(within_cylinder_1, within_cylinder_2)

    Earth_pos_1 = Earth_pos_1[:,within_cylinder]
    Earth_pos_2 = Earth_pos_2[:,within_cylinder]
    
    sol_1 = np.array([x1[within_cylinder], y1[within_cylinder], z1[within_cylinder]])
    sol_2 = np.array([x2[within_cylinder], y2[within_cylinder], z2[within_cylinder]])

#########################################################################################
#Sort cylinder positions based on proximity to starting location on Earth

    dist_1 = np.linalg.norm(Earth_pos_1 - sol_1, axis = 0)
    dist_2 = np.linalg.norm(Earth_pos_1 - sol_2, axis = 0)
    
    Cylinder_pos_1 = np.where(dist_1<dist_2,sol_1, sol_2)
    Cylinder_pos_2 = np.where(dist_1<dist_2,sol_2, sol_1)
    
#########################################################################################
#Recalculate entrance and exit positions on cylinder to account for caps

    Cylinder_pos_1[2][Cylinder_pos_1[2]>=Earth_Radius] = Earth_Radius
    Cylinder_pos_2[2][Cylinder_pos_2[2]>=Earth_Radius] = Earth_Radius
    
    Cylinder_pos_1[2][Cylinder_pos_1[2]<=Earth_Radius-detector_depth] = Earth_Radius-detector_depth
    Cylinder_pos_2[2][Cylinder_pos_2[2]<=Earth_Radius-detector_depth] = Earth_Radius-detector_depth
    
    slope_xz = (Earth_pos_2[2]-Earth_pos_1[2])/(Earth_pos_2[0]-Earth_pos_1[0])
    intercept_xz = Earth_pos_2[2]-slope_xz*Earth_pos_2[0]
    
    Cylinder_pos_1[0], Cylinder_pos_2[0] = (Cylinder_pos_1[2]-intercept_xz)/slope_xz, (Cylinder_pos_2[2]-intercept_xz)/slope_xz
 
    slope_yz = (Earth_pos_2[2]-Earth_pos_1[2])/(Earth_pos_2[1]-Earth_pos_1[1])
    intercept_yz = Earth_pos_2[2]-slope_yz*Earth_pos_2[1]
    
    Cylinder_pos_1[1], Cylinder_pos_2[1] = (Cylinder_pos_1[2]-intercept_yz)/slope_yz, (Cylinder_pos_2[2]-intercept_yz)/slope_yz
    
#########################################################################################
#Filtering trajectories based on an angular cut (most events occur for zenith angles above a given value)
    
    trajectory_direction = Earth_pos_2-Earth_pos_1
    angle = (180/np.pi)*np.arccos(np.einsum('ij,i->j',trajectory_direction, np.array([0,0,1]))/np.linalg.norm(trajectory_direction, axis = 0))
    
    within_angular_range = angle>=angular_cut
    Earth_pos_1 = Earth_pos_1[:,within_angular_range]
    Earth_pos_2 = Earth_pos_2[:,within_angular_range]
    Cylinder_pos_1 = Cylinder_pos_1[:,within_angular_range]
    Cylinder_pos_2 = Cylinder_pos_2[:,within_angular_range]
    
#########################################################################################
#Finding the probability for a neutrino interaction inside the volume (excluding regeneration and NC interactions)
    # trajectory_direction = Earth_pos_2-Earth_pos_1
    # total_L = np.linalg.norm(trajectory_direction, axis = 0)
    # L_ratio_1 = np.linalg.norm(Cylinder_pos_1-Earth_pos_1, axis = 0)
    # L_ratio_2 = np.linalg.norm(Cylinder_pos_2-Earth_pos_1, axis = 0)
    # print(L_ratio_1[0])
    # print(L_ratio_2[0])
    # 
    # print((Cylinder_pos_1-Earth_pos_1)[:,0])
    # print((Cylinder_pos_2-Earth_pos_1)[:,0])
    # 
    # local_normal = Earth_pos_2/np.linalg.norm(Earth_pos_2, axis = 0)
    # theta_EE = (180/np.pi)*np.arccos(np.einsum('ij,ij->j',trajectory_direction, local_normal)/total_L)
    # fig = plt.figure()
    # plt.hist(L_ratio_2-L_ratio_1,bins = 20)
    # 
    # fig = plt.figure()
    # ax = plt.subplot(projection='3d')
    # Xc,Yc,Zc = data_for_cylinder_along_z(0,0,detector_radius,detector_depth, Earth_Radius)
    # ax.plot_surface(Xc, Yc, Zc, color = 'red', alpha=0.5)
    # 
    # # ax.set_xlim3d(-2*detector_radius, 2*detector_radius)
    # # ax.set_ylim3d(-2*detector_radius, 2*detector_radius)
    # # ax.set_zlim3d(Earth_Radius-2*detector_depth, Earth_Radius+2*detector_depth)
    # ax.set_xlim3d(-Earth_Radius, Earth_Radius)
    # ax.set_ylim3d(-Earth_Radius, Earth_Radius)
    # ax.set_zlim3d(-Earth_Radius, Earth_Radius)
    # ax.view_init(elev=0., azim=0)
    # #
    # #for i in range(0,len(Earth_pos_1[0])):
    # for i in range(0, 1):
    #     a = Arrow3D([Earth_pos_1[0][i], Earth_pos_2[0][i]], [Earth_pos_1[1][i], Earth_pos_2[1][i]], [Earth_pos_1[2][i], Earth_pos_2[2][i]], mutation_scale=25, lw=1, arrowstyle="->")
    #     ax.add_artist(a)
    #     #ax.plot([Earth_pos_1[0][i], Earth_pos_2[0][i]], [Earth_pos_1[1][i], Earth_pos_2[1][i]], [Earth_pos_1[2][i], Earth_pos_2[2][i]])
    # #     # ax.scatter([x1[0],x2[0]],[y1[0],y2[0]], [z1[0],z2[0]])
    # #     # ax.scatter(x1[i],y1[i],z1[i], color = 'red')
    # #     # ax.scatter(x2[i],y2[i],z2[i], color = 'green')
    # 
    # ax.scatter(Cylinder_pos_1[0][0], Cylinder_pos_1[1][0], Cylinder_pos_1[2][0], color = 'green')
    # ax.scatter(Cylinder_pos_2[0][0], Cylinder_pos_2[1][0], Cylinder_pos_2[2][0], color = 'red')
    # plt.legend()
    # 
    # plt.show()
    
    # xi, yi = np.meshgrid(energy, angle, indexing='ij')
    # desired_points = np.array([xi.flatten(), yi.flatten()]).T
    # 
    # probabilities = probability_interpolator_function(desired_points)
#########################################################################################
#Append data to the structure array
    
    a1 = np.append(Earth_pos_1, Earth_pos_2, axis = 0)
    b1 = np.append(Cylinder_pos_1, Cylinder_pos_2, axis = 0)
    c1 = np.append(a1, b1, axis = 0) 
    events = np.append(events, c1, axis = 1)

#########################################################################################
#Save data

np.savetxt(str(thrown_events)+"_Example_Trajectories.csv", events.T, delimiter=",")
