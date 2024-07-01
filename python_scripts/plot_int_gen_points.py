import matplotlib.pyplot as plt
import numpy as np


vertex_point=[117842.946250,  -292340.227033,  637583855.087712]
dir_vector=[  0.551554,  0.052592, -0.832480]
test_point=[  -551436122.518741,  -52884413.641951,  1470063462.642013]

enter_point=[ -25313.5, -305991, 6.378e+08 ]
exit_point=[ 5.85709e+08, 5.55453e+07, -2.46269e+08]

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
 
# Creating plot
ax.scatter3D(enter_point[0], enter_point[1],enter_point[2] , color = "green")
ax.scatter3D(exit_point[0], exit_point[1],exit_point[2] , color = "red")
ax.scatter3D(test_point[0], test_point[1],test_point[2] , color = "purple")
ax.scatter3D(vertex_point[0], vertex_point[1],vertex_point[2] , color = "blue")

ax.quiver(vertex_point[0],vertex_point[1],vertex_point[2],dir_vector[0],dir_vector[1],dir_vector[2],length=3e8)
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = 6.378e8*np.cos(u)*np.sin(v)
y = 6.378e8*np.sin(u)*np.sin(v)
z = 6.378e8*np.cos(v)
ax.plot_wireframe(x, y, z, color="r",linewidth=1)
plt.show()