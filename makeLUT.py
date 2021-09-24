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
    num_type=[]
    num_CC = []
    num_NC = []
    num_decays = []
    num_particles = []
    num_gen=[]
    energy = []
    start_energy=[]
    end_pos=[]
    anti_type=[]
    num_GR=[]

    for line in open(filename,mode='r',encoding='utf-8'):
        
        if(lc!=0 and lc!=1 and 'END' not in line):
            num_type.append(int(line.split()[0]))
            anti_type.append(int(line.split()[1]))
            num_NC.append(int(line.split()[2]))
            num_CC.append(int(line.split()[3]))
            num_decays.append(int(line.split()[4]))
            num_GR.append(int(line.split()[6]))
            num_gen.append(int(line.split()[7]))
            #num_particles.append(int(line.split()[5]))
            energy.append(float(line.split()[8]))
            start_energy.append(float(line.split()[9]))
            #end_pos.append(float(line.split()[9]))
            #print line
        lc+=1
    return np.array(num_type),np.array(anti_type),np.array(num_CC), np.array(num_NC), np.array(num_decays), np.array(num_particles), np.array(energy)

tag_list = [tag]
#for thickness in ['0.0', '1.0', '2.0', '3.0', '4.0']:
#  for CS in ['low','mid', 'upp']:
#    #for EL in ['std', 'low']:
#    for EL in ['std', 'low']:
#       tag_list.append('%skm_ice_%sCS_%sEL'%(thickness, CS, EL))
missing_count = 0
#print (tag_list)
for tag in tag_list:
    data_dir = 'testing/particles'
    paths=os.listdir(data_dir)
    eees=[]
    all_energies=['11.0','12.0','13.0','14.0','15.0','16.0','17.0','18.0','19.0','20.0','21.0']
    for e in all_energies:
        for i in paths:
            print(data_dir+'particles_'+e+'_91.4.dat')
            if('particles_'+e+'_91.4.dat'==i):
                eees.append(e)
                break
    
    
    count = 0
    count_true=0
    ang_array = np.concatenate([ np.arange(90.0,95.0,0.1) , np.arange(95.0,181.0,1.0) ])
    th_exit_array = 90.0-ang_array
    e_array = np.array([1e15, 3e15, 1e16, 3e16, 1e17, 3e17, 1e18, 3e18, 1e19, 3e19, 1e20, 3e20, 1e21])
    e_array=np.array([1e15,1e16,1e17,1e18,1e19,1e20,1e21])
   
    outdir = data_dir+'LUT'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    for e in eees:
        #e_s = '%1.0e'%e
        #e_s = e_s.replace('+','')
        #print e_s
        mean_num_CC=[]
        mean_num_NC=[]
        mean_num_decays=[]
        mean_num_gen=[]
        data_array  = []
        type_array=[]
        anti=[]
        #e_pow=log10(e)
        
        for ang in ang_array:
           
            fnm = data_dir+"particles_%s_%.1f.dat"%(e, ang)
            print (fnm)
      
            #num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
            #print(num_CC,num_NC,num_decays,num_particles,energy)
            
            if not os.path.exists(fnm):
                print (fnm, os.path.exists(fnm))
                missing_count += 1
                data_array.append([0])
                type_array.append([0])
                mean_num_CC.append([0])
                mean_num_NC.append([0])
                mean_num_decays.append([0])
                anti.append([0])
            if os.path.exists(fnm):
                the_type,anti_t,num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
                print (fnm, np.size(energy))
                type_array.append(the_type)
                anti.append(anti_t)
                data_array.append(energy)
                mean_num_CC.append(np.mean(num_CC))
                mean_num_NC.append(np.mean(num_NC))
                mean_num_decays.append(np.mean(num_decays))
            np.savez('%s/LUT_%s_eV.npz'%(outdir,e),type_array=type_array, anti=anti,data_array = data_array, th_exit_array = th_exit_array,mean_num_CC=mean_num_CC,mean_num_NC=mean_num_NC,mean_num_decays=mean_num_decays)
            
    
        

print ('missing_count', missing_count)    
"""  
fnm = '/home/ryan/software/NuTauSim/data/Emerging_tau_info_reg_1e+20_90.000000_4.0km_ice_midCS_stdEL.dat'
print('reading file',fnm)
num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
print(num_CC,num_NC,num_decays,num_particles,energy)
"""    
    
