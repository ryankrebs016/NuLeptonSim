#code used to convert NuLeptomSim output to that used for ARA sec. sim

import numpy as np
import sys
part_type=sys.argv[1]
energy=sys.argv[2]
in_file='ARA_finfin/events/'+part_type+'_events_'+energy+'.dat'
out_file='ARA_finfin/with_channel/'+part_type+'_events_'+energy+'.npz'
n_thrown=1000
t_thrown=59068000000
rx=[]
ry=[]
rz=[]
u=[]
v=[]
w=[]
e=[]
y=[]
p_type=[]
i_type=[]
t_num=[]
n_num=[]
regen=[]
channel=[]



with open(in_file,'r') as f:
    line=f.readline()
    while("END" not in line):
        
        arr=line.split(',')
        if(arr[0]=='vert_x'):
            line=f.readline()
            continue
        rx.append(float(arr[0])*10.**-5)
        
       
        ry.append(float(arr[1])*10.**-5)
        rz.append(float(arr[2])*(10.**-5))
        u.append(float(arr[3]))
        v.append(float(arr[4]))
        w.append(float(arr[5]))
        e.append(float(arr[6]))
        y.append(float(arr[7]))
        p_type.append(int(arr[8]))
        i_type.append(int(arr[9]))
        channel.append(int(arr[10]))
        regen.append(int(arr[11]))
        t_num.append(int(arr[12]))
        n_num.append(int(arr[13]))
        line=f.readline()
        
weights=np.ones(np.size(rx))/n_thrown/t_thrown

rx=np.array(rx)
ry=np.array(ry)
rz=np.array(rz)
u=np.array(u)
v=np.array(v)
w=np.array(w)
e=np.array(e)
y=np.array(y)
p_type=np.array(p_type)
i_type=np.array(i_type)
t_num=np.array(t_num)
n_num=np.array(n_num)
channel=np.array(channel)
np.savez(out_file,vert_x=rx,vert_y=ry,vert_z=rz,dir_x=u,dir_y=v,dir_z=w,energy=e,y=y,p_type=p_type,i_type=i_type,weight=weights,regen_code=regen,t_num=t_num,n_num=n_num,channel=channel)


