#code used to convert NuLeptonSim output to that used for ARA sec. sim

import numpy as np
import sys
import os
part_type=sys.argv[1]
energy=sys.argv[2]

#print(part_type,energy)

dir='ARA_FINAL_6/'

in_file=dir+'events/'+part_type+'_events_'+energy+'.dat'
if not os.path.exists(dir+'out_events/'):
    os.mkdir(dir+'out_events/')
out_file=dir+'out_events/'+part_type+'_events_'+energy+'.npz'

#lep_in_file=dir+'leptons/'+part_type+'_leptons_'+energy+'.dat'
#lep_out_file=dir+'out_events/'+part_type+'_leptons_'+energy+'.npz'

t_frac=.1 #n_traj/100000
n_frac=10 #n_throw/1000

n_thrown=1000*n_frac
t_thrown=59068000000*t_frac

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
nc_num=[]
dc_num=[]
gr_num=[]
sto_num=[]

count=0
print('starting events')
with open(in_file,'r') as f:
    line=f.readline()
    while("END" not in line):
        count+=1
        if(count%100000==0): print(count)
        arr=line.split(',')
        if(arr[0]=='vert_x'):
            line=f.readline()
            continue
        try:
            rx.append(float(arr[0])*10.**-5)
            ry.append(float(arr[1])*10.**-5)
            rz.append(float(arr[2])*10.**-5)
            u.append(float(arr[3]))
            v.append(float(arr[4]))
            w.append(float(arr[5]))

            e.append(float(arr[6]))
            y.append(float(arr[7]))

            p_type.append(int(arr[8]))
            i_type.append(int(arr[9]))
            channel.append(int(arr[10]))

            nc_num.append(int(arr[11]))
            dc_num.append(int(arr[12]))
            gr_num.append(int(arr[13]))

            t_num.append(int(arr[14]))
            n_num.append(int(arr[15]))
            sto_num.append(int(arr[16]))

            line=f.readline()
        except:
            print('failed to decode line ',count-1)
            exit()
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
dc_num=np.array(dc_num)
nc_num=np.array(nc_num)
gr_num=np.array(gr_num)
sto_num=np.array(sto_num)
t_num=np.array(t_num)
n_num=np.array(n_num)
channel=np.array(channel)
np.savez(out_file,vert_x=rx,vert_y=ry,vert_z=rz,dir_x=u,dir_y=v,dir_z=w,energy=e,y=y,p_type=p_type,i_type=i_type,weight=weights,nc_num=nc_num,gr_num=gr_num,dc_num=dc_num,sto_num=sto_num,t_num=t_num,n_num=n_num,channel=channel)
#the rest is for the old lepton files in which lep and neutirno events were handfled seperately
print('finished processing')
exit()


print('saved neutrino events')
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
channel=[]
regen=[]
t_num=[]
n_num=[]
ei=[]
ef=[]
frac=[]
sto_type=[]
sto_num=[]
traj=[]
p_num=[]
count=0
print('starting stochastic events')
with open(lep_in_file,'r') as f:
    line=f.readline()
    while ('END' not in line):
        count+=1
        if(count%100000==0): print(count)
        arr=line.split(',')
        if(arr[0]=='logEi'):
            line=f.readline()
            continue
        ei.append(arr[0])
        ef.append(arr[2])
        frac.append(arr[1])
        rx.append(float(arr[3])*10**-5)
        ry.append(float(arr[4])*10**-5)
        rz.append(float(arr[5])*10**-5)
        u.append(arr[6])
        v.append(arr[7])
        w.append(arr[8])
        sto_type.append(arr[9])
        traj.append(arr[10])
        p_num.append(arr[11])
        sto_num.append(arr[12])
        line=f.readline()
sto_text=[]
weights=np.ones(np.size(ei))/n_thrown/t_thrown
for i in range(len(sto_type)):
    if(sto_type[i]==0):
        sto_text.append('brem')
    if(sto_type[i]==1):
        sto_text.append('pp')
    if(sto_type[i]==2):
        sto_text.append('pn')
            	
np.savez(lep_out_file,weights=weights,vert_x=rx,vert_y=ry,vert_z=rz,dir_x=u,dir_y=v,dir_z=w,ei=ei,edep=ef,e_frac=frac,sto_type=sto_text,t_num=traj,n_num=p_num,sto_num=sto_num)
print('ending stochastic events')
