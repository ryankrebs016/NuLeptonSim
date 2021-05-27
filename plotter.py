#plotter plots energies. plots plots exit probs
import matplotlib.pyplot as plt
import os
import numpy as np

LUTdir='icecube/'

energy_list=['1e+15','1e+16','1e+17','1e+18','1e+19']
part_index=[1,2,4,5]
names=["nuMu","nuTau","Muon","Tau"]


def process_lut_for_parts(LUT_fnm, particle_type): #same as load LUT but looks for specific particles. IE want just taus

    f = np.load(LUT_fnm,allow_pickle=True)
    type_array=f['type_array']
    th_array = f['th_exit_array']
    th_em_array = -th_array
    data_array = f['data_array']
    P_exit = np.zeros(len(th_array))
    mean_num_CC     = f['mean_num_CC']
    mean_num_NC     = f['mean_num_NC']
    mean_num_decays = f['mean_num_decays']
    proc_type=[]
    proc_energy=[]
    #print(np.size(th_array),np.size(data_array),np.size(type_array))
   

    for k in range(0,len(th_array)):
        #print(np.size(type_array[k]),np.size(data_array[k]))
        temp_energy=[]
        temp_type=[]
        if(len(data_array[k])!=0 ):
            part_count=0
            #print(np.size(type_array[k]))
            #print(np.size(data_array[k]))
            for i in range(0,len(data_array[k])):
                #print(type_array[k][i],data_array[k][i])
                #type_array[k][i]
                if(type_array[k][i]==particle_type):
                    temp_type.append(type_array[k][i])
                    temp_energy.append(data_array[k][i])
                    part_count+=1
            proc_type.append(temp_type)
            proc_energy.append(temp_energy)
               
            

            P_exit[k] = part_count/1e4

    return type_array,th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays
#returned elements shouldn't have any other particle than type 2 - nutau

#energy='1e+16'
def plot_per_particle():
    for i in part_index:
        plt.figure(1,figsize=(12,8))
        plt.xscale('log')
        print(i)
        plt.ylabel('exit energy log(E/eV)')
        plt.xlabel('exit angle (deg)')
        for energy in ['1e+16','1e+17','1e+18','1e+19']:

            print(LUTdir+'LUT_%s_eV.npz'%energy)
            mean_array=[]
            std_array=[]
            ang_array=[]
            ang_num=0
            if os.path.exists(LUTdir+'LUT_%s_eV.npz'%(energy)) and os.path.getsize(LUTdir+'LUT_%s_eV.npz'%(energy)) > 0:
                type_array,th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = process_lut_for_parts(LUTdir+'/LUT_%s_eV.npz'%(energy),i)
                for j in range(np.size(th_em_array)):
                    #mean_energy=np.mean(th_em_array[])
                    if(P_exit[j]>0):
                        mean_array.append(np.mean(data_array[j]))
                        std_array.append(np.std(data_array[j]))
                        ang_array.append(th_em_array[j])           
            plt.errorbar(ang_array,mean_array,std_array,capsize=0.2,elinewidth=0.2)
        plt.axis([1,100,14,20])       
        plt.legend(['1e+16','1e+17','1e+18','1e+19'])    
        plt.savefig(LUTdir+'energy_dist_of_type_%s'%i)
        plt.close()
                
def plot_per_energy():
    for energy in ['1e+16','1e+17','1e+18','1e+19']:
        print(energy)
        plt.figure(1,figsize=(12,8))
        plt.xscale('log')
        plt.grid(True)
        
        plt.ylabel('exit energy log(E/eV)')
        plt.xlabel('exit angle (deg)')
        for i in part_index:
            print(i)
            print(LUTdir+'LUT_%s_eV.npz'%energy)
            mean_array=[]
            std_array=[]
            ang_array=[]
            ang_num=0
            if os.path.exists(LUTdir+'LUT_%s_eV.npz'%(energy)) and os.path.getsize(LUTdir+'LUT_%s_eV.npz'%(energy)) > 0:
                type_array,th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = process_lut_for_parts(LUTdir+'/LUT_%s_eV.npz'%(energy),i)
                for j in range(np.size(th_em_array)):
                    #mean_energy=np.mean(th_em_array[])
                    if(P_exit[j]>0):
                        mean_array.append(np.mean(data_array[j]))
                        std_array.append(np.std(data_array[j]))
                        ang_array.append(th_em_array[j])           
            plt.errorbar(ang_array,mean_array,std_array,capsize=0.2,elinewidth=0.2)
        plt.axis([1,100,14,20])         
        plt.legend(part_index)    
        plt.savefig(LUTdir+'particle_dist_of_energy_%s'%energy)
        plt.close()



def main():
    plot_per_energy()       
    plot_per_particle()        



if __name__ == "__main__":

    main()


