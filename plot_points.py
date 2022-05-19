import numpy as np
import matplotlib.pyplot as plt
import pandas


xi=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[0],dtype=float))
yi=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[1],dtype=float))
zi=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[2],dtype=float))
xf=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[3],dtype=float))
yf=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[4],dtype=float))
zf=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[5],dtype=float))
xe=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[6],dtype=float))
ye=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[7],dtype=float))
ze=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[8],dtype=float))

xex=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[9],dtype=float))
yex=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[10],dtype=float))
zex=np.array(pandas.read_csv('in_files/59068000000_Example_Trajectories.csv',usecols=[11],dtype=float))

path_length=[]
zeniths=[]
x_vec=[]
y_vec=[]
z_vec=[]
for i in range(np.size(xi)):
    len=((xex[i]-xi[i])**2+(yex[i]-yi[i])**2+(zex[i]-zi[i])**2)**(1/2)
    x_vec.append((xex[i]-xi[i])/len)
    y_vec.append((yex[i]-yi[i])/len)
    z_vec.append((zex[i]-zi[i])/len)
    path_length.append(len)
    zeniths.append(180/np.pi*np.arccos((zf[i]-zi[i])/len))
    if i%10000==0:print(i)
filtered=[]
#for i in range(np.size(zeniths)):
#    if zeniths[i]<90:
#        filtered.append(zeniths[i])

plt.figure(2)
plt.hist(np.array(x_vec))
plt.yscale("log")
plt.ylim([1000,100000])
plt.title("x")

plt.figure(3)
plt.hist(np.array(y_vec))
plt.yscale("log")
plt.title("y")
plt.ylim([1000,100000])
plt.figure(4)
plt.hist(np.array(z_vec))
plt.yscale("log")
plt.title("z")
plt.ylim([1000,100000])
plt.figure(5)
plt.hist(np.array(zeniths))
plt.yscale("log")
plt.title("zen")
plt.show()
exit()
plt.figure(3)
plt.hist(np.array(path_length))
plt.show()
fig1=plt.figure(1)
ax1=fig1.add_subplot(projection='3d')
ax1.scatter3D(xex,yex,zex,color='red')

#ax.scatter3D(xi,yi,zi)
fig2=plt.figure(2)
ax2=fig2.add_subplot(projection='3d')
ax2.scatter3D(xe,ye,ze,color='green')
plt.show()
