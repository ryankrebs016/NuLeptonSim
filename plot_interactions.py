#makes plots of histograms from event generator

from lzma import CHECK_NONE
from re import X
from tkinter import Y
import numpy as np
import pandas
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
doubles=0
energy='21'
x=[]
y=[]
z=[]
u=[]
v=[]
w=[]
int_type=[]
zeniths=[]
throws=[]
traj=[]
db_traj=[]
db_count=0
ey=[]
e=[]
event_num=0
rads=[]
with open('ARA_finfin/events/nutau_events_%s.0.dat'%energy,'r') as f:
    line=f.readline()
    temp_traj=0
    while("END" not in line):
        
        arr=line.split(',')
        if(arr[0]=='vert_x'):
            line=f.readline()
            continue
        x.append(float(arr[0])*10**-5)
        y.append(float(arr[1])*10**-5)
        z.append(float(arr[2])*(10**-5)-6378)
        u.append(float(arr[3])*1)
        v.append(float(arr[4])*1)
        w.append(float(arr[5])*1)
        e.append(float(arr[6]))
        ey.append(float(arr[7]))
        rads.append(((x[event_num])**2+(y[event_num])**2)**(1/2))
        zeniths.append(180/np.pi*np.arccos(float(arr[5])))
        int_type.append(int(arr[9]))
        throws.append(int(arr[11]))
        traj.append(int(arr[10]))
        if int(arr[10]) == temp_traj: doubles=doubles+1
        temp_traj=int(arr[10])
        event_num=event_num+1
        line=f.readline()
for i in range(np.size(throws)-1):
    if(throws[i-1]==throws[i] and throws[i-1]==throws[i]):
        #print("bang bang")
        db_count+=1
print(doubles," multiple events in same traj")
print(event_num," total events")
print('total db ',db_count)
print('total events',event_num)
zeniths=np.array(zeniths)
color_int=[]
int_type_name=[]
temp_zen=[]
nu_zen=[]
tau_zen=[]
weighted_events=[0,0,0,0,0]

bin=[80]
for i in range(20):
    bin.append(i*5+85)
    bin.append(i*5+85)

#print(bin)
total_cc=0
cc_e=[]
nc_e=[]
tau_x=[]
tau_y=[]
tau_z=[]
tau_u=[]
tau_v=[]
tau_w=[]
tau_e=[]
tau_ey=[]
tau_rad=[]
cc_y=[]
nc_y=[]
count_zeros=0

for i in range(np.size(int_type)):
    if(int_type[i]==0):
        color_int.append('b')
        int_type_name.append('CC')
        weighted_events[int_type[i]]=weighted_events[int_type[i]]+1/1000
        nu_zen.append(zeniths[i])
        total_cc+=1
        if np.log10(e[i])+9 == 21:
            cc_e.append(np.log10(e[i])+9)
            cc_y.append(ey[i])
    elif(int_type[i]==1):
        color_int.append('r')
        int_type_name.append('NC')
        weighted_events[int_type[i]]=weighted_events[int_type[i]]+1/1000
        nu_zen.append(zeniths[i])
        nc_e.append(np.log10(e[i])+9)
        nc_y.append(ey[i])
    elif(int_type[i]==2):
        color_int.append('m')
        int_type_name.append('GR_h')
        temp_zen.append(zeniths[i])
        weighted_events[int_type[i]]=weighted_events[int_type[i]]+1/1000
    elif(int_type[i]==3):
        color_int.append('k')
        int_type_name.append('GR_l')
        weighted_events[int_type[i]]=weighted_events[int_type[i]]+1/1000
    elif(int_type[i]==4):
        color_int.append('g')
        int_type_name.append('D')
        weighted_events[int_type[i]]=weighted_events[int_type[i]]+1/1000
        tau_zen.append(zeniths[i])
        tau_e.append(np.log10(e[i])+9)
        tau_ey.append(ey[i])
        tau_x.append(x[i])
        tau_y.append(y[i])
        tau_z.append(z[i])
        tau_u.append(u[i])
        tau_v.append(v[i])
        tau_w.append(w[i])
        tau_rad.append(rads[i])
        if ey[i]==0:
            count_zeros=count_zeros+1
    else:
        color_int.append('y')

'''
plt.figure(1)
plt.hist(np.array(x),histtype='step',density=True)
plt.hist(np.array(tau_x),histtype='step',density=True)
plt.ylabel("count")
plt.xlabel("x")
plt.yscale('log')
plt.legend(["all","taus"])

plt.figure(2)
plt.hist(np.array(y),histtype='step',density=True)
plt.hist(np.array(tau_y),histtype='step',density=True)
plt.ylabel("count")
plt.xlabel("y")
plt.yscale('log')
plt.legend(["all","taus"])

plt.figure(7)
plt.hist(np.array(rads),histtype='step',density=True)
plt.hist(np.array(tau_rad),histtype='step',density=True)
plt.ylabel("count")
plt.xlabel("rad")
plt.yscale('log')
plt.legend(["all","taus"])
plt.figure(3)
plt.hist(np.array(z),histtype='step',density=True)
plt.hist(np.array(tau_z),histtype='step',density=True,bins=20)
plt.ylabel("count")
plt.xlabel("z")
plt.yscale('log')
plt.legend(["all","taus"])
plt.figure(4)
plt.hist(np.array(u),histtype='step',density=True)
plt.hist(np.array(tau_u),histtype='step',density=True)
plt.ylabel("count")
plt.xlabel("u")
plt.yscale('log')
plt.legend(["all","taus"])
plt.figure(5)
plt.hist(np.array(v),histtype='step',density=True)
plt.hist(np.array(tau_v),histtype='step',density=True)
plt.ylabel("count")
plt.xlabel("v")
plt.yscale('log')
plt.legend(["all","taus"])
plt.figure(6)
plt.hist(np.array(w),histtype='step',density=True)
plt.hist(np.array(tau_w),histtype='step',density=True,bins=20)
plt.ylabel("count")
plt.xlabel("w")
plt.yscale('log')
plt.legend(["all","taus"])
'''
plt.figure(8)
plt.hist2d(w,rads,bins=20,density=True)
plt.xlabel("w")
plt.ylabel("radius")
plt.title("all %s"%energy)
plt.colorbar()
plt.figure(9)
plt.hist2d(tau_w,tau_rad,bins=20,density=True)
plt.xlabel("w")
plt.ylabel("radius")
plt.title("tau %s"%energy)
plt.colorbar()
plt.show()


print(len(nu_zen),len(tau_zen))
exit()

plt.figure()
plt.hist(tau_ey,density=True)
print(len(tau_ey))
print(count_zeros)
plt.show()

exit()
   
plt.figure(4)
plt.bar(['CC','NC','GR_h','GR_l','D'],weighted_events)
plt.yscale('log')
plt.ylim([.01,10000])
plt.title("10^%s eV weighted events"%energy)
weight=np.ones(np.size(np.array(zeniths)))*1/1000
#plt.figure(2)
#plt.hist(np.array(int_type_name))
#plt.yscale('log')
#plt.title('Event type histogram')
plt.savefig('ARA_fin/plots/%s_events.png'%energy)
plt.figure(3)
plt.hist(np.array(nu_zen),weights=weight[0:np.size(np.array(nu_zen))],bins=bin,histtype='step')
plt.hist(np.array(tau_zen),weights=weight[0:np.size(np.array(tau_zen))],bins=bin,histtype='step')
plt.legend(['neutrinos events','tau events'])
plt.xlim([80,180])
plt.xlabel("zenith angle")
plt.title('10^%s eV nu event trajectory zeniths'%energy)
plt.show()
#plt.savefig('ARA_fin/plots/%s_zeniths.png'%energy)
exit()
plt.figure(5)
plt.hist(np.array(tau_zen),weights=weight[0:np.size(np.array(tau_zen))],bins=bin)
plt.xlim([80,180])
plt.xlabel("zenith angle")
plt.title('10^17 eV tau event trajectory zeniths')

print(np.mean(np.array(zeniths)))
#plt.savefig('ARA/plots/21_zen.png')
plt.show()
exit()
fig=plt.figure(1)
ax=fig.add_subplot(111, projection='3d')

ax.quiver(x[0:500],y[0:500],z[0:500],u[0:500],v[0:500],w[0:500])
ax.set_xlim3d(-15,15)
ax.set_ylim3d(-15,15)
ax.set_zlim3d(-15,15)
plt.show()
