from contextlib import nullcontext
import numpy as np
import matplotlib as plt
from numpy.lib.npyio import save
import scipy as sc
import scipy.integrate
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit
import pandas as pd
import time
import pickle
from matplotlib.animation import FuncAnimation

process=True

dir='testing_GR/'
particle_dir=dir+'particles/'
event_dir=dir+'events/'
bin_dir=dir+'binned/'
LUTdir=dir+'LUT/'
prob_dir=dir+'p_exit/'
obj_dir='obj/'
#set bin energies to something that works for all energy ranges. Easier to interpolate the missing energies

def main():
    if(len(sys.argv)>1):
        energy_list=[sys.argv[1]]
    else: energy_list=['11.0','12.0','13.0','14.0','15.0','16.0','17.0','18.0','19.0','20.0','21.0']
    prep_files(6)
    prep_files(5)
    prep_files(4)
    #calc_eff_taus()
    #calc_neutrino_eff_area()
    #get_E_from_CC(18)
    #calc_eff_area()
    #plot_dN_muon(19.0)

    
    #plot_an_angle('19.0')
    #
    # plot_an_angle('18.0')
    """
    for i in range(10):
        plot_dN_tau(20.0,i*10+1,'muon')

    
    plot_an_angle('21.0')
    """
    exit()
    #plot_p_exit([11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0],'nutau')
    #exit()
    

    #calc_eff_area()
    '''
    #prep_files(5) #5 for taus, 4 for muons
    plot_p_exit(['14.0'],'nutau')
    plot_p_exit(['16.0'],'nutau')
    plot_p_exit(['18.0'],'nutau')
    plot_p_exit(['20.0'],'nutau')
    plt.legend(['14','16','18','20'])
    plt.xlabel('emergence angle')
    plt.ylabel('tau exit prob')
    plt.show()
    '''
    
    
    angles=5,10,15,20,25,30,36,40
    
    #for i in angles:
        #plot_dN_tau('18.0',i,'tau')
    
    
    to_plot=1
    
    if(os.path.exists(dir+'effective_areas.npz') and to_plot==True):

        f=np.load('final_muon/effective_areas.npz',allow_pickle=True)
        x=f['energies']
        y=f['effective_areas']
        p=np.load('final_tau/effective_areas.npz',allow_pickle=True)
        a=p['energies']
        b=p['effective_areas']
        
        r=np.load('final_tau/nu_effective_areas.npz',allow_pickle=True)
        t=r['energies']
        u=r['effective_areas']
        ccs=[]
        loccs=[]
        '''
         for i in range(np.size(u)):
            u[i]=u[i]*cross_section(10**t[i])*5*10**28*1/54
            ccs.append(cross_section(10**t[i]))
            loccs.append(np.log10(cross_section(10**t[i])))
            print(u[i])
            print(ccs)
            print(loccs)
     
        
        '''
       
      
        ff=np.load('no_regen_tau/effective_areas.npz',allow_pickle=True)
        tt=ff['energies']
        uu=ff['effective_areas']
        

        ice_energies=np.array(pd.read_csv('icecube_tau_exposure.csv',usecols=[0]))
        ice_exposures=np.array(pd.read_csv('icecube_tau_exposure.csv',usecols=[1]))
       
        auger_energies=np.array(pd.read_csv('icecube_auger_sensitivity.csv',usecols=[0]))
        auger_exposures=np.array(pd.read_csv('icecube_auger_sensitivity.csv',usecols=[1]))
       
        ic80_allsky_t_andm_E=np.array(pd.read_csv('IC80_allsky_tandm.csv',usecols=[0]))
        ic80_allsky_t_andm_eff=np.array(pd.read_csv('IC80_allsky_tandm.csv',usecols=[1]))

        ic40_ps_allsky_E=np.array(pd.read_csv('IC40_PS_allsky.csv',usecols=[0]))
        ic40_ps_allsky_eff=np.array(pd.read_csv('IC40_PS_allsky.csv',usecols=[1]))

        ic80_and_deepcore_E=np.array(pd.read_csv('IC80_and_deepcore.csv',usecols=[0]))
        ic80_and_deepcore_eff=np.array(pd.read_csv('IC80_and_deepcore.csv',usecols=[1]))

        for i in range(np.size(ice_energies)):
            ice_energies[i]=ice_energies[i]+9
            ice_exposures[i]=10**ice_exposures[i]/((100*100)*2426*24*60*60*4*np.pi)
        for i in range(np.size(auger_energies)):
            auger_energies[i]=auger_energies[i]+9
            auger_exposures[i]=10**auger_exposures[i]/((100*100)*2426*24*60*60*4*np.pi)

        for i in range(np.size(ic40_ps_allsky_E)):
            ic40_ps_allsky_E[i]=ic40_ps_allsky_E[i]+9
            ic40_ps_allsky_eff[i]=10** ic40_ps_allsky_eff[i]
        
        for i in range(np.size(ic80_allsky_t_andm_E)):
            ic80_allsky_t_andm_E[i]=ic80_allsky_t_andm_E[i]+9
            ic80_allsky_t_andm_eff[i]=10**ic80_allsky_t_andm_eff[i]
        
        for i in range(np.size(ic80_and_deepcore_E)):
            ic80_and_deepcore_E[i]=ic80_and_deepcore_E[i]+9
            ic80_and_deepcore_eff[i]=10**ic80_and_deepcore_eff[i]
        total=[]
        for i in range(np.size(t)):
            total.append(y[i]+b[i]+u[i])
        plt.figure(1)
        plt.title('Icecube Effective Area')
        plt.xlabel('log(E_v/eV)')
        plt.ylabel('effective area m^2')
        #plt.semilogy(x,y,linewidth=3)
        plt.semilogy(a,b,linewidth=3)
        #plt.semilogy(tt,uu,linewidth=3)
        #plt.semilogy(t,u,linewidth=3)
        #plt.semilogy(t,total)
        #plt.semilogy(ice_energies,ice_exposures)
        plt.semilogy(auger_energies,auger_exposures)
        #plt.semilogy(ic40_ps_allsky_E,ic40_ps_allsky_eff)
        #plt.semilogy(ic80_allsky_t_andm_E,ic80_allsky_t_andm_eff)
        #plt.semilogy(ic80_and_deepcore_E,ic80_and_deepcore_eff)
        plt.xlim(12,21)
        plt.ylim(1,10**5)
        plt.legend(['NLS Vt+Vm','IC 80 all sky Vt+Vm'])
        #plt.legend(['NuLeptonSim Vm','NuLeptonSim Vt Decay','NuLeptonSim Vt cascade','Icecube Vt (2018)','IC 80 All Sky Vt and Vm'])
        plt.grid(True,which='both')
        plt.show()
        exit()
    
    #prep_files()
    
    '''
    for i in range(50):
        aeff=[]
        for j in range(50):
            print(i)
            aeff.append(calc_effective_area(19,t_energies[i],m_energies[j],1))
        aeffs.append(aeff)
    
    plt.figure(1)
    for i in range(50):
        plt.plot(t_energies,aeffs[i])
    
    plt.show()
    '''
    #plot_an_angle('19.0')
def plot_an_angle(energy,what_plot='counts'):
    fig=plt.figure()
    ax=plt.axes(xlim=(13,22),ylim=(0,5))
    ax.set_xlabel('Log(E_t / eV)')
    ax.set_ylabel('dN/dE')
    ax.set_title('NuTau energy ditributions from 10^19 eV neutrino\nFrom 0 degrees to 90 degrees')
    line,=ax.plot([],[],lw=2)
    print('plotting')
    angle_text=ax.text(.02,.95,'',transform=ax.transAxes)
    number_text=ax.text(.02,.80,'',transform=ax.transAxes)
    th_em_array,bins_by_ang,derivs_by_ang,middle_bins=read_Ntau_files(bin_dir+str(energy)+'_binned_nutau.npz')
    def init():
        line.set_data([],[])
        angle_text.set_text('')
        number_text.set_text('')
        return line,

    def animate(i):
       
        norm,mid_bins=get_normalized(bins_by_ang,middle_bins,i,'tau')
        x=mid_bins
        y=norm
        line.set_data(x,y)
        angle_text.set_text('angle %.1f degrees'%th_em_array[i])
        number_text.set_text('%s particles'%np.size(norm))
        return line,

    anim=FuncAnimation(fig,animate,init_func=init,frames=135,interval=1)
    anim.save('test_plots/19nuTau_energy_distribution.gif', fps=30)
    plt.show()
    
    # get rid from here 
    
    x=np.linspace(11,21,np.size(bins_by_ang[0]))
    plt.figure(1,figsize=[10,8])
    plt.title('10^'+energy+' eV Neutrino Energy. '+'nuEl energy distribution')
    plt.xlabel("log(E)")

    if(what_plot=='counts'):
        plt.ylabel("normalized dN/dE")
        for i in range(10):
            norm,mid_bins=get_normalized(bins_by_ang,middle_bins,i,'tau')
            plt.plot(mid_bins,norm)


    plt.legend(th_em_array[25::10])
    plt.show()

    
def cross_section(E_v):
    c0=-1.826
    c1=-17.31
    c2=-6.406
    c3=1.431
    c4=-17.91
    eps=np.log10(E_v*10**-9)
    print(eps)
    logcc=c1+c2*np.log(eps-c0)+c3*(np.log(eps-c0))**2+c4/np.log(eps-c0)
    cc=10**logcc
    return cc
class get_E_from_CC():
    
    def __init__(self):
        self.store={'null':[0]}
        if(os.path.exists(obj_dir+'nutau_decays.pkl')):
            self.load_energies()
         
    def load_energies(self):
        with open(obj_dir+'nutau_decays.pkl','rb') as f:
            return pickle.load(f)
    def save_energies(self):
        with open(obj_dir+'nutau_decays.pkl','wb') as f:
            pickle.dump(self.store,f,pickle.HIGHEST_PROTOCOL)
    def check_storage(self,E,E_t):
        E_G=E-9  
        self.jline=0
        self.jmin=int(10000*(E_G-1)*1)
        self.jmax=int(self.jmin+1001) 
        e_min=E-6
        e_max=E
        bin_num=50
        e_slope=(e_max-e_min)/bin_num
        bin_energies=np.linspace(e_min,e_max,bin_num)
        try:
            norm_bins=self.store['%s'%E]
            _,dN=general_find_and_interpolate(bin_energies,norm_bins,E_t,interp=True)
            return dN
        except:
            e_in=np.array(pd.read_csv('tables/final_cteq5_cc_nu.data',header=1,usecols=[1],delimiter=' ',dtype='float'))
            e=np.array(pd.read_csv('tables/final_cteq5_cc_nu.data',header=1,usecols=[2],delimiter=' ',dtype='float'))
            Ei=[]
            #print(e)
            #print(np.size(e))
            
            Eo=e[self.jmin+1:self.jmax:1]
            '''
            for i in range(np.size(e_in)):
                if(self.jmin<i and i<self.jmax):
                    #print(e_in[i])
                    Eo.append(e[i])
            print(np.size(Eo),np.size(Etest))
            print(Etest[0],Etest[np.size(Etest)-1])
            print(Eo[0],Eo[np.size(Eo)-1])
            '''
            e_min=E-6
            e_max=E
            bin_num=100
            e_slope=(e_max-e_min)/bin_num
            bin_energies=np.linspace(e_min,e_max,bin_num)
            bins=np.zeros(bin_num)
            for i in range(np.size(Eo)):
                Eo[i]=np.log10(Eo[i]*10**E)
                index=int((Eo[i]-e_min)/e_slope)
                if(index<0):
                    index=0
                if(index>bin_num):
                    index=bin_num-1
                bins[index]=bins[index]+1
            norm_bins,_=get_normalized(bins,bin_energies,0,'nu')
            self.store.update({'%s'%E:norm_bins})
            plt.figure(1)
            plt.plot(bin_energies,bins)
            plt.show()
            plt.figure(1)
            plt.plot(bin_energies,norm_bins)
            plt.show()
            print(sc.integrate.trapz(norm_bins,bin_energies))

            _,dN=general_find_and_interpolate(bin_energies,norm_bins,E_t,interp=True)
            return dN
        


    '''

    
    
    r = []
    fy = open('tables/final_cteq5_cc_nu.data')
    
    
 
    
    e_in=pd.read_csv('tables/final_cteq5_cc_nu.data',header=0,usecols=[1],delimiter=' ',dtype='float')
    e=pd.read_csv('tables/final_cteq5_cc_nu.data',header=0,usecols=[2],delimiter=' ',dtype='float')
    print(e_in,e)
    exit()
	
    for line in fy:
        print(line)
        if line.strip():      # Read only non-empty lines, skip empty lines 
            yy = line.split() # Split lines using empty space as separator
	
            if yy[0]!='#':    #  only read lines not starting with '#'
                jline = jline + 1
	
                if jmin < jline < jmax:  # read only lines between jmin and jmax
                    yvalue = float(yy[1])  # 2nd column (index 1) is the Bjorken y
                    r.append(yvalue)

    fy.close()
   
    dy = 0.02
    bins_y = np.arange(0., 1.0+dy, dy)
    #print(E*r[10])
    for i in range(np.size(r)):
        r[i]=np.log10(r[i]*10**E)
        index=int((r[i]-e_min)/e_slope)
        if(index<0):
            index=0
        if(index>bin_num):
            index=bin_num-1
        bins[index]=bins[index]+1
    '''
  

#general use interpolation function. takes just scalars
def general_interp_value(x,y,i):
    slope=(y[1]-y[0])/(x[1]-x[0])
    return slope*(i-x[0])+y[0]

#general use function that finds first index for interpolation and then interpolates to return index and the number. takes 1d vectors
def general_find_and_interpolate(x,y,i,interp=False):

    for j in range(np.size(x)):
        x[j]=float(x[j])
    index=0
    value=y[0]
    out_of_domain=False
    
    if(i<=np.min(x)):
        #print('i below range')
        out_of_domain=True
        if interp==True:
            return index,y[index]
        return index, value
    if(i>=np.max(x)):
        #print('i above range')
        index=np.size(x)-1
        out_of_domain=True
        if interp==True:
            return index,y[index-1]
        return index, value
    #print('i in range')
    for j in range(np.size(x)-1):
        
        if(i >= x[j] and i < x[j+1]):
            index=j
            break
            
    
    if(interp==True and out_of_domain==False):
        value=general_interp_value([x[index],x[index+1]],[y[index],y[index+1]],i)
    
    return index,value

#defines the bin energies for individual bins from the min, max, bin size, and index
def def_bin_energy(min_energy,max_energy,bin_size,index):

    bin_min=index*bin_size+min_energy
    bin_max=(index+1)*bin_size+min_energy

    return bin_min,bin_max

#performs a finite diffeence derivitive. takes 1d vectors for x and y
def discrete_der(middle_bins,bins):

    derivs=[]
    for i in range(np.size(bins)-1):
        derivs.append((bins[i+1]-bins[i])/(middle_bins[i+1]-middle_bins[i]))
   
    return derivs

#opens and gets info from dN/dE tau files
def read_Ntau_files(filename):

    f=np.load(filename,allow_pickle=True)
    th_em_array=f['th_em_array']
    bins_by_ang=f['bins_by_ang']
    derivs_by_ang=f['derivs_by_ang']
    middle_bins=f['middle_bins']

    return th_em_array,bins_by_ang,derivs_by_ang,middle_bins

#gets info then calculates the normalized dN/dE tau and interpolates a value
def dN_by_angle(v_energy,t_energy,angle,type):
    
    filename=bin_dir+str(v_energy)+'_binned_%s.npz'%type
    th_em_array,bins_by_ang,derivs_by_ang,middle_bins=read_Ntau_files(filename)
    #testing to see if this makes it better norm_dervs,mid_bins=get_normalized(derivs_by_ang,middle_bins,angle,'tau')
    norm,mid_bins=get_normalized(bins_by_ang,middle_bins,angle,type)
    
    which_dN,dN=general_find_and_interpolate(mid_bins,norm,t_energy,True)
    
    return float(dN)

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
        nutau_count=0
        tau_count=0
        muon_count=0
        numu_count=0
        nuel_count=0
        if(len(data_array[k])!=0 ):
            part_count=0
            
            #print(np.size(type_array[k]))
            #print(np.size(data_array[k]))
            for i in range(0,len(data_array[k])):
                #print(type_array[k][i],data_array[k][i])
                #type_array[k][i]
                if(type_array[k][i]==1):
                    nuel_count+=1
                if(type_array[k][i]==2):
                    numu_count+=1
                if(type_array[k][i]==3):
                    nutau_count+=1
                if(type_array[k][i]==5):
                    muon_count+=1
                if(type_array[k][i]==6):
                    tau_count+=1

                if(type_array[k][i]==particle_type):
                    temp_type.append(type_array[k][i])
                    temp_energy.append(data_array[k][i])
                    part_count+=1
        
            proc_type.append(temp_type)
            proc_energy.append(temp_energy)
              
            

            P_exit[k] = part_count/1e6
            print(part_count)
        #print(np.size(temp_type))     
        #print(numu_count,muon_count,nutau_count,tau_count)
    return proc_type,th_em_array, P_exit, proc_energy, mean_num_CC, mean_num_NC, mean_num_decays
#returned elements shouldn't have any other particle than the type specified
    

def bin_energies(energy,part,bin_num,LUTdir):
   
    type_array,th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = process_lut_for_parts(LUTdir+'LUT_%s_eV.npz'%(energy),part)
    #print(LUTdir+'LUT_%s_eV.npz'%(energy))
    bins_by_ang=[]
    derivs_by_ang=[]
    for ang_ind in range(np.size(th_em_array)):
        try:
            print(np.shape(data_array),' ',ang_ind)
            energies=data_array[ang_ind]
            max_energy=22. #max_energy=np.max(energies)
            min_energy=10. #min_energy=np.min(energies)
            bin_size=(max_energy-min_energy)/bin_num
            bins=np.zeros(bin_num)
            min_bin=np.zeros(bin_num)
            max_bin=np.zeros(bin_num)
            middle_bins=[]
            
            for j in range(bin_num):
                bin_min,bin_max=def_bin_energy(min_energy,max_energy,bin_size,j)
                #print(j, "    ",bin_min,"    ", bin_max)
                min_bin[j]=bin_min
                max_bin[j]=bin_max
                #print(np.size(energies))
                middle_bins.append((max_bin[j]-min_bin[j])/2+min_bin[j])
            for i in range(np.size(energies)):
                index=int(bin_num*(energies[i]-min_energy)/(max_energy-min_energy))
                if(index<0):
                    print('particle below min energy')
                    index=0
            
                if(index>=bin_num):
                    print('particle above max energy',energies[i],th_em_array[ang_ind],energy,part)
                    
                    
                    index=bin_num-1
                bins[index]=bins[index]+1
                #print('index is %s'%index)
                
            derivs=middle_bins #discrete_der(middle_bins,bins)       essentially forget aboutdoing d derivs 
            bins_by_ang.append(bins)
            derivs_by_ang.append(derivs)
        except:
            bins_by_ang.append(np.zeros(bin_num))
            derivs_by_ang.append(np.zeros(bin_num))
    save_string=''
    if(part==6):
        save_string='tau'
    if(part==3):
        save_string='nutau'
    if(part==5):
        save_string='muon'
    if(part==2):
        save_string='numu'
    if(part==1):
        save_string='nuel'
    if(part==4):
        save_string='el'
    print(bin_dir+str(energy)+"_binned_%s.npz"%save_string)
    np.savez(bin_dir+str(energy)+"_binned_%s.npz"%save_string,bin_low=min_bin,bin_high=max_bin,th_em_array=th_em_array,bins_by_ang=bins_by_ang,derivs_by_ang=bins_by_ang, middle_bins=middle_bins)

def interpolate_counts(energy,bin_dir):
    energies=np.array([12,13,14,15,16,17,18,19])
    interpolate=True
    energy_index_range=0
    for i in range(np.size(energies)-1):
        if(energy<energies[i+1] and energy>=energies[i]):
            energy_index_range=i
        if(energy==energies[i]):
            interpolate=False
    th_em_array1,bins_by_ang1,derivs_by_ang1=read_files(bin_dir+energy_index_range)
    if(interpolate==True):
        th_em_array2,bins_by_ang2,derivs_by_ang2=read_files(bin_dir+energy_index_range+1)
        interpolated_bins=[]
        interpolated_derivs=[]
        what_value=(energy-energies[energy_index_range])/(energies[energy_index_range+1]-energies[energy_index_range])*100
        #force what value to .01 accuracy since i will interpolate 100 points between the energies
        for i in range(np.size(th_em_array1)):
            xvals=np.linspace(energies[energy_index_range],energies[energy_index_range],100)
            counts=np.interpolate(xvals,[energies[energy_index_range],energies[energy_index_range+1]],[bins_by_ang1[i],bins_by_ang2[i]])
            dervs=np.interpolate(xvals,[energies[energy_index_range],energies[energy_index_range+1]],[derivs_by_ang1[i],derivs_by_ang2[i]])
            interpolated_bins.append(counts)
            interpolated_derivs.append(dervs)

        this_one=0        
        for i in range(100):
            if what_value>=i and what_value==i+1:
                this_one=i
        final_counts=interpolated_bins[::1][this_one]
        final_derivs=interpolated_derivs[::1][this_one]
        return th_em_array1,final_counts,final_derivs
    print('no interpolation')
    return th_em_array1,bins_by_ang1,derivs_by_ang1
    
def get_muon_dist():
    f=open('pythia_tau.txt','r')
    rows=[]
    for i in range(100000):
        rows.append(f.readline())
   
    rel_energies=[]
    for i in range(np.size(rows)):
        if(rows[i].find(':13,')!=-1):
            rel_energy=float(rows[i][rows[i].find(':13,')+4:rows[i].find(':-14,')])
            #print(rel_energy)
            rel_energies.append(rel_energy)

    #print(np.size(rel_energies))
    #print(np.min(np.array(rel_energies)))
    #print(np.max(np.array(rel_energies)))
    np.savez('muon_energies.npz',rel_energies=rel_energies)
    #print(rel_energies)

def bin_muon_energies(filename):
    f=np.load(filename,allow_pickle=True)
    rel_energies=f['rel_energies']
    print(rel_energies)
    for i in range(np.size(rel_energies)):
        rel_energies[i]=np.log10(rel_energies[i])
    print(rel_energies)
    
    min_energy=np.log10(.0001)
    print(min_energy)
    
    max_energy=np.log10(1.01)
    print(max_energy)
    bin_num=100
    bin_size=(max_energy-min_energy)/bin_num
    bins=np.zeros(bin_num)
    min_bin=np.zeros(bin_num)
    max_bin=np.zeros(bin_num)
    middle_bins=np.zeros(bin_num)
    for j in range(bin_num):
        bin_min,bin_max=def_bin_energy(min_energy,max_energy,bin_size,j)
            #print(j, "    ",bin_min,"    ", bin_max)
        min_bin[j]=bin_min
        max_bin[j]=bin_max
        middle_bins[j]=(bin_max-bin_min)/2+bin_min
            #print(np.size(energies))

        for i in range(np.size(rel_energies)):
            if(rel_energies[i]>bin_min and rel_energies[i]<bin_max):
                bins[j]=bins[j]+1
   
    
    np.savez('muon_dist.npz',middle_bins=middle_bins,bins=bins)
    
 
    
def get_muon_dN(t_energy,m_energy):
    f=np.load('muon_dist.npz')
    middle_bins=f['middle_bins']
    bins=f['bins']
    #print(np.size(middle_bins),np.size(bins))
    
    xax=np.zeros(np.size(middle_bins))
    for i in range(np.size(xax)):
        xax[i]=t_energy+middle_bins[i]
    derivs=discrete_der(xax,bins)
    norm_bins,_=get_normalized(bins,xax,12,'mu')
    index,dN=general_find_and_interpolate(xax,norm_bins,m_energy,True)
    return dN

def interp_aeff(filename,muon_energy):
    f=open(filename)
    energies=[]
    aeffs=[]
    for i in range(18):
        row=f.readline()
        energies.append(9+float(row[:row.find(',')]))
        aeffs.append(10**float(row[row.find(',')+1:row.find('\n')]))
    '''
    plt.figure(1)
    plt.title('icecube muon effective area')
    plt.xlabel('muon energy log(E/GeV)')
    plt.ylabel('effective area m^2')
    plt.plot(energies,aeffs)
    plt.show()
    '''
    
    if muon_energy < energies[0]:
        #print(' energy below')
        return aeffs[0]
    if muon_energy> energies[np.size(energies)-1]:
        #print('energy above')
        return aeffs[np.size(aeffs)-1]
    
    index,aeff=general_find_and_interpolate(energies,aeffs,muon_energy,True)
 
    return aeff
    
#normalizes deriv vs energy distributions
def get_normalized(derivs, middle_bins,angle,type):
    #print(middle_bins,derivs)

    if(type=='tau' or type=='nutau' or type=='muon' or type=='numu'):
        #index by angle
        dervs=derivs[int(angle)]
        middles=middle_bins
   
    if(type=='mu' or type=='nu'):
        #dont index by angle
        dervs=derivs
        middles=middle_bins
    dervs_size=np.size(dervs)
    middles_size=np.size(middles)
    #print(dervs_size,middles_size)
   
    while dervs_size!=middles_size:
        if(dervs_size>middles_size):
            dervs=dervs[:dervs_size-1:1]
        if(dervs_size<middles_size):
            middles=middles[:middles_size-1:1]   
        dervs_size=np.size(dervs)
        middles_size=np.size(middles)
    #print(dervs,middles)
    area=sc.integrate.trapz(dervs,middles)
 
  
    #print('area is %s'%area)
    normalized_dervs=np.zeros(np.size(dervs))
    if(area!=0):
        for i in range(np.size(dervs)):
            normalized_dervs[i]=dervs[i]/area

  
    return normalized_dervs,middle_bins

#processes the LUT to make a smaller file for just probabilites
def process_P(v_energy,type='tau'):
    #print(v_energy)
    part_num=6
    other='nutau'
    if(type=='muon'):
        part_num=5
        other='numu'
    if(type=='electron'):
        part_num=4
        other='nuel'
    nu_num=part_num-3
    print(prob_dir+v_energy+'_p_exit_%s.npz'%type)
    #print(LUTdir+'LUT_%s_eV.npz'%(v_energy))

    type_array,th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = process_lut_for_parts(LUTdir+'LUT_%s_eV.npz'%(v_energy),part_num)
    np.savez(prob_dir+v_energy+'_p_exit_%s.npz'%type,th_em_array=th_em_array,P_exit=P_exit)

    type_array,th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = process_lut_for_parts(LUTdir+'LUT_%s_eV.npz'%(v_energy),nu_num)
    np.savez(prob_dir+v_energy+'_p_exit_%s.npz'%other,th_em_array=th_em_array,P_exit=P_exit)


#returns the individual number for the probability 
def find_p(filename,angle):
    f=np.load(filename)
    P_exit=f['P_exit']
    th_em_array=f['th_em_array']
    index,prob=general_find_and_interpolate(th_em_array,P_exit,angle,True)
    return prob
    
#input variables in terms of log(E/eV)... will return the effective area (integrand in the big integral for the decay path)
def calc_effective_area(v_energy,t_energy,m_energy,angle): #used for just tau decay processes
    print('calculating tau eff area for v energy: %s, angle: %s, t energy: %s, and mu energy: %s'%(v_energy,angle,t_energy,m_energy))
    sys.stdout.write("\033[F")

    #find Psurv
    #print('finding p surv')
    th_em_array,_,_,_=read_Ntau_files(bin_dir+'14.0_binned_tau.npz')
    filename_energy_list=[11,12,13,14,15,16,17,18,19,20,21]
    filename_energies=['1e+11','1e+12','1e+13','1e+14','1e+15','1e+16','1e+17''1e+18','1e+19','1e+20','1e21']
    log_energies=['11','12','13','14','15','16','17','18','19','20','21']
    which_energy, _ = general_find_and_interpolate(log_energies,[0],v_energy)

    p1=find_p(prob_dir+'%s_p_exit_tau.npz'%log_energies[which_energy],angle)
    prob=p1
    if(which_energy+1<np.size(log_energies)):

        p2=find_p(prob_dir+'%s_p_exit_tau.npz'%log_energies[which_energy+1],angle)
        prob=general_interp_value([log_energies[which_energy],log_energies[which_energy+1]],[p1,p2],v_energy)
  
    if(prob==0):return 0
    #print('p surv is %s\n'%prob)

    #find dN/dE tau
    #print('finding dN tau')
    which_v,tmp=general_find_and_interpolate(log_energies,[0],v_energy)
    which_angle,tmp=general_find_and_interpolate(th_em_array,[0],angle)
    use_energy=log_energies[which_v]
    use_angle=th_em_array[angle]
    which_angle=1
    dN11=dN_by_angle(log_energies[which_v],t_energy,angle,'tau') #th_em_array[which_angle] already interpolates on t_energy so just need to do it for v energy
    dN_tau=dN11
    if(which_energy+1<np.size(log_energies)):
        dN12=dN_by_angle(log_energies[which_v+1],t_energy,angle,'tau')
        dN_tau=general_interp_value([log_energies[which_energy],log_energies[which_energy+1]],[dN11,dN12],v_energy)
    
    if dN_tau==0:return 0
    #print('dN tau is %s\n'%dN11)

    #find dN/dE muon
    #print('getting dN mu')
   
    dN_mu=get_muon_dN(t_energy,m_energy)
    #print('dN_mu is %s\n'%dN_mu)
    if(dN_mu==0): return 0

    #set branching ratio
    branch=.18
  
    #find aeff of muon events
    #print('getting aeff of icecube for muon events')
    aeff=interp_aeff('Default Dataset.csv',m_energy)
    #print('aeff of icecube is %s\n'%aeff)

    effective_area=prob*dN_tau*dN_mu*branch*aeff
    #print('effective area is %s\n\n'%effective_area)
    return effective_area

    
def calc_deff_area_n(CC,v_energy,v_prime_energy,t_energy,m_energy,angle):

    Np=1 #rough estimation using area of proton and a 1km^3 icecube volume
    th_em_array,_,_,_=read_Ntau_files(bin_dir+'14.0_binned_tau.npz')
    filename_energy_list=[11,12,13,14,15,16,17,18,19,20,21]
    filename_energies=['1e+11','1e+12','1e+13','1e+14','1e+15','1e+16','1e+17''1e+18','1e+19','1e+20','1e21']
    log_energies=['11','12','13','14','15','16','17','18','19','20','21']
    which_energy, _ = general_find_and_interpolate(log_energies,[0],v_energy)

    #calc p
    p1=find_p(prob_dir+'%s_p_exit_nutau.npz'%log_energies[which_energy],angle)
    prob=p1
    if(which_energy+1<np.size(log_energies)):
        p2=find_p(prob_dir+'%s_p_exit_nutau.npz'%log_energies[which_energy+1],angle)
        prob=general_interp_value([log_energies[which_energy],log_energies[which_energy+1]],[p1,p2],v_energy)
    if(prob==0): 
        return 0

    #calc dN_v
    dN11=dN_by_angle(log_energies[which_energy],v_prime_energy,angle,'numu') #th_em_array[which_angle] already interpolates on t_energy so just need to do it for v energy
    dN_v=dN11
    if(which_energy+1<np.size(log_energies)):
        dN12=dN_by_angle(log_energies[which_energy+1],v_prime_energy,angle,'numu')
        dN_v=general_interp_value([log_energies[which_energy],log_energies[which_energy+1]],[dN11,dN12],v_prime_energy)
    
    if dN_v==0:
        return 0

    
    #calc dN_m
    dN31=get_muon_dN(t_energy,m_energy)
    if(dN31==0): 
        return 0

    #calc dN_t
    dN21=CC.check_storage(v_prime_energy,t_energy)
    if(dN21==0): 
        return 0
    #get icecube detector level muon effecive area
    aeff=interp_aeff('Default Dataset.csv',m_energy)

    effective_area=prob*dN_v*Np*dN21*dN31*.18*aeff
    return effective_area

def calc_neutrino_eff_area():
    v_energies=np.arange(11,21.5,1)
    v_prime_energies=np.arange(11,21.5,.2)
    t_energies=np.arange(11,21.5,.2)
    m_energies=np.arange(11,21.5,.2)
    angles_indicies=np.arange(0,135,1,dtype='int')
    ind_vars={'v_energies':v_energies,'v_prime_energies':v_prime_energies,'ang_ind':angles_indicies,'m_energies':m_energies,'t_energies':t_energies}
    if(os.path.exists(dir+'raw_eff_areas.npz')and False):
        f=np.load(dir+'raw_eff_areas.npz',allow_pickle=True)
        daeff=f['daeff']
        print(daeff)
        integrate_neutrino_stuff(ind_vars,daeff)
    daeff=[]
    CC=get_E_from_CC()
    for v in v_energies:
        aeff_v=[]
        #print(v)
        for ang in angles_indicies:
            aeff_ang=[]
            for v_prime in v_prime_energies:
                if(v_prime>v):
                    break
                aeff_v_prime=[]
                for t in t_energies:
                    if(t>v_prime):
                        break
                    aeff_t=[]
                    for m in m_energies:
                        if(m>t):
                            break
                        print('nu tau eff area for v energy: %s, ang ind: %s, v prime energy: %s, t energy: %s, m energy: %s'%(v,ang,v_prime,t,m))
                        sys.stdout.write("\033[F")
                        aeff_m=calc_deff_area_n(CC,v,v_prime,t,m,ang)
                        aeff_t.append(aeff_m)
                    aeff_v_prime.append(aeff_t)
                aeff_ang.append(aeff_v_prime)
            aeff_v.append(aeff_ang)
        daeff.append(aeff_v)
    CC.save_energies()
    try:
        integrate_neutrino_stuff(ind_vars,daeff)
    except: 
        np.savez(dir+'raw_eff_areas.npz',daeff=daeff,ind_vars=ind_vars,dtype='object')
        exit()

def integrate_neutrino_stuff(ind_vars,daeff):
    th_em_array,_,_,_=read_Ntau_files(bin_dir+'14.0_binned_tau.npz')
    aeff=[]
    for v in range(np.size(ind_vars['v_energies'])):
        for theta in range(np.size(daeff[v])):
            for v_prime in range(np.size(daeff[v][theta])):
                for t in range(np.size(daeff[v][theta][v_prime])):
                    #print(v,theta,v_prime,t)
                    #print(daeff[v][theta][v_prime][t])
                    if(np.size(daeff[v][theta][v_prime][t])==1):
                        area1=0
                    else:
                        area1=sc.integrate.trapz(daeff[v][theta][v_prime][t],x=ind_vars['m_energies'][0:np.size(daeff[v][theta][v_prime][t]):1])
                    #print(area1)
                    
                    
                    daeff[v][theta][v_prime][t]=area1
                area2=sc.integrate.trapz(daeff[v][theta][v_prime],x=ind_vars['t_energies'][0:np.size(daeff[v][theta][v_prime]):1])
                daeff[v][theta][v_prime]=area2
            area3=sc.integrate.trapz(daeff[v][theta],x=ind_vars['v_prime_energies'][0:np.size(daeff[v][theta]):1])
            daeff[v][theta]=area3
        area4=sc.integrate.trapz(daeff[v],x=th_em_array[0:np.size(daeff[v]):1])
        aeff.append(area4/2)
    diff=np.size(ind_vars['v_energies'])-np.size(aeff)
    adjust_x=0
    adjust_y=0
    if diff>0:
        adjust_x=diff
        adjust_y=0
        print('diff non zero')
    if diff<0:
        adjust_x=0
        adjust_y=-diff
        print('diff non zero')
    np.savez(dir+'nu_effective_areas.npz',energies=ind_vars['v_energies'][0:np.size(ind_vars['v_energies'])-adjust_x:1],effective_areas=aeff[0:np.size(aeff)-adjust_y:1])
    
def plot_p_exit(energy,type):
    plt.figure(1)
    for e in energy:
        f=np.load(prob_dir+'%s_p_exit_%s.npz'%(e,type))
        P_exit=f['P_exit']
        th_em_array=f['th_em_array']
        plt.loglog(th_em_array,P_exit)
    plt.ylim(10**-7,1)
    plt.legend(energy)
    plt.grid(True)
    plt.ylabel('Exit Probability')
    plt.xlabel('Emergence Angle (degrees)')
    plt.show()

def plot_dN_tau(v_energy,angle,type):
    t_energy=np.linspace(11,21,100)
    dN=np.zeros(100)
    sums=[]
    for i in range(100):
        dN[i]=(dN_by_angle(v_energy,t_energy[i],angle,type))
    plt.figure(1)
    plt.plot(t_energy,dN)
    plt.xlabel('log(E_tau/eV)')
    plt.ylabel('normalized dN/dE')
    plt.title('normalized dN/dE for %s'%type)
    plt.show()

def plot_dN_muon(t_energy):
    m_energies=np.linspace(11,21,50)
    dN=[]
    for i in range(50):
        dN.append(get_muon_dN(t_energy,m_energies[i]))
    plt.figure(1)
    plt.plot(m_energies,dN)
    plt.xlabel('log(E_muon/eV)')
    plt.ylabel('normalized dN/dE')
    plt.title('normalized dN/dE for muons')
    plt.show()


def test_interpolation():
    a=[3,4,5]
    b=[7,8,9]
    x=2
    y=5
    z=6
    print('low test')
    print(general_find_and_interpolate(a,b,x,True))    
    print('mid test')
    print(general_find_and_interpolate(a,b,y,True))    
    print('high test')
    print(general_find_and_interpolate(a,b,z,True))    
      
'''
def aeff_objective(a,b,c,d,x):
    return a*(b*x+c)**(1/2)+d

def extrap_aeff():

    f=open('Default Dataset.csv')
    energies=[]
    aeffs=[]
    for i in range(18):
        row=f.readline()
        energies.append(float(row[:row.find(',')]))
        aeffs.append(10**float(row[row.find(',')+1:row.find('\n')]))
    print(aeffs)
    popt,_=curve_fit(aeff_objective,energies,aeffs)  
    a,b,c,d=popt
    print('%s*(%s*x+%s)**(1/2)+%s'%(a,b,c,d)) 

def aeff_func(energy):
    y= 0.06574498*(0.018130432*energy+5.1445770)**(1/2)+0.065744983
    y=176603*(118*energy+55852)**(1/2)+189257
    return y
def extrapolate_aeff_muon(energies):
    y=[]
    for i in range(np.size(energies)):
        y.append(i)
    return y
'''
def prep_files(type):
    if(not os.path.exists(bin_dir)):
        os.mkdir(bin_dir)
    if(not os.path.exists(prob_dir)):
        os.mkdir(prob_dir)    
    if(not os.path.exists('muon_energies.npz')):
        get_muon_dist()
        bin_muon_energies('muon_energies.npz')

    energy_list=['11.0','12.0','13.0','14.0','15.0','16.0','17.0','18.0','19.0','20.0','21.0','22.0']
    energy_list=['15.6','15.7','15.9']

    part_tau=type #taus
    part_nu=type-3
    bin_num=1000
    if(type==6):
        type_str='tau'
    if(type==5):
        type_str='muon'
    if(type==4):
        type_str='electron'
    for energy in energy_list:
        #print(LUTdir+'LUT_%s_eV.npz'%(energy))
        if(os.path.exists(LUTdir+'LUT_%s_eV.npz'%(energy))):
            print('processing probs')
            process_P(energy,type_str)
            print('processing taus')
            _=bin_energies(energy,part_tau,bin_num,LUTdir)
            print('processing nutaus')
            _=bin_energies(energy,part_nu,bin_num,LUTdir)

def plot_effective_area_by_angle(aeff):
    
    angles=[]
    
    #print('aeff is ',aeff)
    #print(np.shape(aeff))
    #print(np.size(aeff))
    #print('aeff[0] is ',aeff[0])
    #print(aeff[1])
    plt.figure(1)
    th_em_array,_,_,_=read_Ntau_files(bin_dir+'14.0_binned_muon.npz')
    for i in range(np.size(aeff,axis=0)):
        plt.semilogy(th_em_array,aeff[i])
    plt.legend(np.arange(11,21.5,1))
    plt.xlabel('emergence angle')
    plt.ylabel('Aeff(E_v,theta)')
    plt.show()
    plt.figure(1)
    th_em_array,_,_,_=read_Ntau_files(bin_dir+'14.0_binned_muon.npz')
    for i in range(np.size(aeff,axis=0)):
        plt.loglog(th_em_array,aeff[i])
    plt.legend(np.arange(11,21.5,1))
    plt.xlabel('emergence angle')
    plt.ylabel('Aeff(E_v,theta)')
    plt.show()



def integrate_stuff(e_ind,th_ind,ind_vars,d_aeff):
    aeff=[]
    th_em_array,_,_,_=read_Ntau_files(bin_dir+'14.0_binned_tau.npz')
    #print(e_ind,th_ind,ind_vars)
    aeff_by_angle=[]
    plotting=[]
    for v in range(np.size(e_ind)):
        for theta in range(np.size(th_ind)):
            for t in range(np.size(d_aeff[v][theta])):
                area1=sc.integrate.trapz(d_aeff[v][theta][t],x=ind_vars['m_energies'][0:np.size(d_aeff[v][theta][t]):1])
        
                if(area1<0): 
                    print('area1 is neg',area1)
                    print(np.sum(np.array(d_aeff[v][theta][t])))
                    print(np.array(ind_vars['m_energies'][0:np.size(d_aeff[v][theta][t]):1]))
                    print(np.array(d_aeff[v][theta][t]))
                    exit()
                d_aeff[v][theta][t]=area1
                
            area2=sc.integrate.trapz(d_aeff[v][theta],x=ind_vars['t_energies'][0:np.size(d_aeff[v][theta]):1])
            
            if(area2<0): 
                print('area2 is neg')
                exit()
            d_aeff[v][theta]=area2
        plotting.append(d_aeff[v])
        
        area3=2*np.pi*sc.integrate.trapz(d_aeff[v],x=th_em_array[0:np.size(d_aeff[v]):1])
        print(th_em_array[0:np.size(d_aeff[v]):1])
        if(area3<0): 
            print('area3 is neg')
            exit()
        d_aeff[v]=area3/(4*np.pi)
    #print(d_aeff)
    #plot_effective_area_by_angle(plotting)
    diff=np.size(ind_vars['v_energies'])-np.size(d_aeff)
    adjust_x=0
    adjust_y=0
    if diff>0:
        adjust_x=diff
        adjust_y=0
    if diff<0:
        adjust_x=0
        adjust_y=-diff

    np.savez(dir+'effective_areas.npz',energies=ind_vars['v_energies'][0:np.size(ind_vars['v_energies'])-adjust_x:1],effective_areas=d_aeff[0:np.size(d_aeff)-adjust_y:1])
    '''
    plt.figure(1)
    
    
    
    plt.semilogy(ind_vars['v_energies'][0:np.size(ind_vars['v_energies'])-adjust_x:1],d_aeff[0:np.size(d_aeff)-adjust_y:1])
    plt.show()
    '''
    
def calc_deff_area_of_muon(v_energy,theta_index,m_energy): #used for only muons 
    log_energies=['11','12','13','14','15','16','17','18','19','20','21']
    which_energy, _ = general_find_and_interpolate(log_energies,[0],v_energy)
    #find p exit
    #print('finding p')

    p1=find_p(prob_dir+'%s_p_exit_muon.npz'%log_energies[which_energy],theta_index)
    prob=p1
    
    if(which_energy+1<np.size(log_energies)):

        p2=find_p(prob_dir+'%s_p_exit_muon.npz'%log_energies[which_energy+1],theta_index)
        prob=general_interp_value([log_energies[which_energy],log_energies[which_energy+1]],[p1,p2],v_energy)
  
    if(prob==0):return 0

    #find dN_mu
    #print('finding dN')
    dN11=dN_by_angle(log_energies[which_energy],m_energy,theta_index,'muon') #th_em_array[which_angle] already interpolates on t_energy so just need to do it for v energy
    dN_muon=dN11

    if(which_energy+1<np.size(log_energies)):
        dN12=dN_by_angle(log_energies[which_energy+1],m_energy,theta_index,'muon')
        dN_muon=general_interp_value([log_energies[which_energy],log_energies[which_energy+1]],[dN11,dN12],v_energy)

    if dN_muon==0:return 0
    if(v_energy>16.3 and v_energy<16.7):
        print('P1 is %s, P2 is %s, dN1 is %s, dN2 is %s'%(p1,p2,dN11,dN12))
        print('P is %s, dN is %s'%(prob,dN_muon))
        exit()
    #find Aeff mu
    #print('find aeff icecube')
    aeff=interp_aeff('Default Dataset.csv',m_energy)

    effective_area=prob*dN_muon*aeff
    #print(prob,dN_muon,aeff,effective_area)
 
    return effective_area

def calc_eff_area():
    time_start=time.perf_counter()
    v_energies=np.arange(11,21.5,1)
    m_energies=np.arange(11,21.5,.01)
    #print(v_energies)
    #print(m_energies)
    angles_indicies=np.arange(0,135,1,dtype='int')
    e_ind=np.arange(0,np.size(v_energies),1,dtype='int')
    ind_vars={'v_energies':v_energies,'ang_ind':angles_indicies,'m_energies':m_energies}
    daeff=[]
    for v in v_energies:
        aeff_v=[]
        
        #print(v)
        
        for ang in angles_indicies:
            #print(ang)
            aeff_theta=[]
            for m in m_energies:
                print('v energy: %s, ang ind: %s, m_energy: %s'%(v,ang,m))
                sys.stdout.write("\033[F")
                if(m>v): break
                aeff_m=calc_deff_area_of_muon(v,ang,m)
                aeff_theta.append(aeff_m)
            aeff_v.append(aeff_theta)
        daeff.append(aeff_v)
    integrate_muon_stuff(e_ind,angles_indicies,ind_vars,daeff)
    elapsed_time=(time.perf_counter()-time_start)/60
    print(elapsed_time)

def integrate_muon_stuff(e_ind,angles_indicies,ind_vars,daeff):
    th_em_array,_,_,_=read_Ntau_files(bin_dir+'14.0_binned_muon.npz')
    aeff=[]
    #print(th_em_array)
   
    plotting=[]
    for v in range(np.size(e_ind)):
        i=0
        
        for theta in range(np.size(angles_indicies)):
            
            area1=sc.integrate.trapz(daeff[v][theta],x=ind_vars['m_energies'][0:np.size(daeff[v][theta]):1])
            
            
            #if(area1>0):
                #plt.figure(1)
                #print('area 1 is ',area1)
                #plt.plot(ind_vars['m_energies'][0:np.size(daeff[v][theta]):1],daeff[v][theta][::]) 
                #plt.show()
            daeff[v][theta]=area1
        
        plotting.append(daeff[v])
        area2=2*np.pi*sc.integrate.trapz(daeff[v],x=th_em_array[0:np.size(daeff[v]):1])
        daeff[v]=area2/(4*np.pi)
    plot_effective_area_by_angle(plotting)
    aeff=daeff
    diff=np.size(ind_vars['v_energies'])-np.size(aeff)
    adjust_x=0
    adjust_y=0
    if diff>0:
        adjust_x=diff
        adjust_y=0
        print('diff non zero')
    if diff<0:
        adjust_x=0
        adjust_y=-diff
        print('diff non zero')
    np.savez(dir+'effective_areas.npz',energies=ind_vars['v_energies'][0:np.size(ind_vars['v_energies'])-adjust_x:1],effective_areas=aeff[0:np.size(aeff)-adjust_y:1])

def calc_eff_taus():
    v_energies=np.arange(11,21.5,1)
    thetas=np.arange(0,135,1,dtype='int')
    m_energies=np.arange(11,21.5,.1)
    t_energies=np.arange(11,21.5,.1)
    ind_vars={'v_energies':v_energies,'thetas':thetas,'m_energies':m_energies,'t_energies':t_energies}
    #calc_effective_area(20,19,15,12)
    e_ind=np.arange(0,np.size(v_energies),1,dtype='int')
    th_ind=np.arange(0,np.size(thetas),1,dtype='int')
   
    #plot_dN_muon(19)
    
    
    aeffs=[]
    count=0
    aeff_v=[]
    for v in v_energies:
        aeff_theta=[]
        for theta in thetas:
            aeff_t=[]
            for t in t_energies:
                if t>v: 
                    #print('t>v')
                    break
                aeff_m=[]
                for m in m_energies:
                    #print(v,'\t',theta,'\t',t,'\t',v)
                    

                    if(m>t): 
                        #print('m>t')
                        break
                    eff=calc_effective_area(v,t,m,theta)
                    #if(eff>0):
                        #print(eff)
                        
                    aeff_m.append(eff)
                aeff_t.append(aeff_m)
            aeff_theta.append(aeff_t)
        aeff_v.append(aeff_theta)
    #indexed like aeff_v[v_index][theta_index][t_index][m_index]=effective area for v,theta,t,m
    integrate_stuff(e_ind,th_ind,ind_vars,aeff_v)


if __name__ == "__main__":
    main()

