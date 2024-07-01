import numpy as np

par_dir='ARA_FINAL_2/'
event_dir='out_events/'
broken_dir='broken_up/'
filename_wo='numu_events_21.0'
filename='numu_events_21.0.npz'


data=np.load(par_dir+event_dir+filename)

total_events=np.size(data[data.files[0]])
print(total_events)

for i in range(int(total_events/1E6)):
    s1=int(1E6*i)
    s1=int(int(total_events/1E6)*1E6)
    s2=int(1E6*(i+1))
    s2=total_events
    i=67
    print(i,s1,s2,total_events)
    if i==int(total_events/1E6) and False:
        np.savez(par_dir+broken_dir+filename_wo+'_%i.npz'%i,vert_x=data['vert_x'][s1:total_events:1],vert_y=data['vert_y'][s1:total_events:1],vert_z=data['vert_z'][s1:total_events:1],dir_x=data['dir_x'][s1:total_events:1],dir_y=data['dir_y'][s1:total_events:1],dir_z=data['dir_z'][s1:total_events:1],energy=data['energy'][s1:total_events:1],y=data['y'][s1:total_events:1],p_type=data['p_type'][s1:total_events:1],i_type=data['i_type'][s1:total_events:1],nc_num=data['nc_num'][s1:total_events:1],gr_num=data['gr_num'][s1:total_events:1],dc_num=data['dc_num'][s1:total_events:1],sto_num=data['sto_num'][s1:total_events:1],t_num=data['t_num'][s1:total_events:1],n_num=data['n_num'][s1:total_events:1],channel=data['channel'][s1:total_events:1])
        break
    np.savez(par_dir+broken_dir+filename_wo+'_%i.npz'%i,vert_x=data['vert_x'][s1:s2:1],vert_y=data['vert_y'][s1:s2:1],vert_z=data['vert_z'][s1:s2:1],dir_x=data['dir_x'][s1:s2:1],dir_y=data['dir_y'][s1:s2:1],dir_z=data['dir_z'][s1:s2:1],energy=data['energy'][s1:s2:1],y=data['y'][s1:s2:1],p_type=data['p_type'][s1:s2:1],i_type=data['i_type'][s1:s2:1],nc_num=data['nc_num'][s1:s2:1],gr_num=data['gr_num'][s1:s2:1],dc_num=data['dc_num'][s1:s2:1],sto_num=data['sto_num'][s1:s2:1],t_num=data['t_num'][s1:s2:1],n_num=data['n_num'][s1:s2:1],channel=data['channel'][s1:s2:1])
    exit()

print('done')
