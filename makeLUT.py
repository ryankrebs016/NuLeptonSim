import matplotlib
#matplotlib.use('Agg')
from pylab import *
import os

#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161115/1.5km_ice'
#tag = '1.5km_ice'
#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161115/3.7km_ocn'
#tag = '3.7km_ocn'
#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161128/1.5km_ice_lowCS'
#tag = '1.5km_ice_lowCS'
#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161128/3.7km_ocn_lowCS'
#tag = '3.7km_ocn_lowCS'

tag = '0.0km_ice_midCS_stdEL'
#tag = '0.0km_ice_regCS_lowEL'
#tag = '0.0km_ice_lowCS_regEL'
#tag = '0.0km_ice_lowCS_lowEL'
#tag = '1.0km_ice_regCS_regEL'
#tag = '1.0km_ice_regCS_lowEL'
#tag = '1.0km_ice_lowCS_regEL'
#tag = '1.0km_ice_lowCS_lowEL'
#tag = '2.0km_ice_regCS_regEL'
#tag = '2.0km_ice_regCS_lowEL'
#tag = '2.0km_ice_lowCS_regEL'
#tag = '2.0km_ice_lowCS_lowEL'
#tag = '3.0km_ice_regCS_regEL'
#tag = '3.0km_ice_regCS_lowEL'
#tag = '3.0km_ice_lowCS_regEL'
#tag = '3.0km_ice_lowCS_lowEL'
#tag = '4.0km_ice_regCS_regEL'
#tag = '4.0km_ice_regCS_lowEL'
#tag = '4.0km_ice_lowCS_regEL'
#tag = '4.0km_ice_lowCS_lowEL'

def read_emerging(filename):
    lc = 0
    num_CC = []
    num_NC = []
    num_decays = []
    num_particles = []
    energy = []
    for line in open(filename,mode='r',encoding='utf-8'):
        if(lc!=0 and 'END' not in line):
            num_CC.append(int(line.split()[0]))
            num_NC.append(int(line.split()[1]))
            num_decays.append(int(line.split()[2]))
            num_particles.append(int(line.split()[3]))
            energy.append(float(line.split()[4]))
            #print line
        lc+=1
    return np.array(num_CC), np.array(num_NC), np.array(num_decays), np.array(num_particles), np.array(energy)

tag_list = []
for thickness in ['0.0', '1.0', '2.0', '3.0', '4.0']:
  for CS in ['low','mid', 'upp']:
    #for EL in ['std', 'low']:
    for EL in ['std', 'low']:
        tag_list.append('%skm_ice_%sCS_%sEL'%(thickness, CS, EL))
missing_count = 0
#print (tag_list)
for tag in tag_list:
    data_dir = '/home/ryan/software/NuTauSim/data'
    count = 0
    count_true=0
    ang_array = np.concatenate([ np.arange(90.0,95.0,0.1) , np.arange(95.0,180.0,1.0) ])
    th_exit_array = 90.0-ang_array
    e_array = np.array([1e15, 3e15, 1e16, 3e16, 1e17, 3e17, 1e18, 3e18, 1e19, 3e19, 1e20, 3e20, 1e21])

   
    outdir = '/home/ryan/software/NuTauSim/LUTs/%s'%tag
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    for e in e_array:
        e_s = '%1.0e'%e
        #e_s = e_s.replace('+','')
        #print e_s
        mean_num_CC=[]
        mean_num_NC=[]
        mean_num_decays=[]
        data_array  = []
        
        for ang in ang_array:
           
            fnm = data_dir+'/Emerging_tau_info_reg_%s_%1.1f00000_%s.dat'%(e_s, ang, tag)
            print (fnm)
            #num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
            #print(num_CC,num_NC,num_decays,num_particles,energy)
            
            if not os.path.exists(fnm):
                print (fnm, os.path.exists(fnm))
                missing_count += 1
                data_array.append([0])
            if os.path.exists(fnm):
                num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
                print (fnm, np.size(energy))
                data_array.append(energy)
                mean_num_CC.append(np.mean(num_CC))
                mean_num_NC.append(np.mean(num_NC))
                mean_num_decays.append(np.mean(num_decays))
                np.savez('%s/LUT_%s_eV.npz'%(outdir,e_s), data_array = data_array, th_exit_array = th_exit_array,mean_num_CC=mean_num_CC,mean_num_NC=mean_num_NC,mean_num_decays=mean_num_decays)
            
    
        

print ('missing_count', missing_count)    
"""  
fnm = '/home/ryan/software/NuTauSim/data/Emerging_tau_info_reg_1e+20_90.000000_4.0km_ice_midCS_stdEL.dat'
print('reading file',fnm)
num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
print(num_CC,num_NC,num_decays,num_particles,energy)
"""    
    
