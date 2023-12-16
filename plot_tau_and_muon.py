import numpy as np
import matplotlib.pyplot as plt

energies=np.arange(6,12.5,1)


data_dir="data/"
cont_dist_filename="cont_dist_"
sto_dist_filename="sto_dist_"

muons_c=[]
muons_s=[]

taus_c=[]
taus_s=[]

avg_m_c=[]
avg_m_s=[]
avg_t_c=[]
avg_t_s=[]


for e in energies:
  #muons_c=np.loadtxt(data_dir+cont_dist_filename+e+".0_"+13+".txt")
  #muons_s=np.loadtxt(data_dir+sto_dist_filename+e+".0_"+13+".txt")
  taus_c=np.loadtxt(data_dir+cont_dist_filename+str(e)+".0_"+str(15)+".txt")[0]*10**-5
  taus_s=np.loadtxt(data_dir+sto_dist_filename+str(e)+".0_"+str(15)+".txt")[0]*10**-5

  avg_t_c.append([np.mean(taus_c),np.std(taus_c)])
  avg_t_s.append([np.mean(taus_s),np.std(taus_s)])
  
  #avg_m_c.append([np.mean(muons_c),np.std(muons_c)])
  #avg_m_s.append([np.mean(muons_s),np.std(muons_s)])


plt.figure()
plt.errorbar(energies,avg_t_c[:][0],yerr=avg_t_c[:][1])
plt.errorbar(energies,avg_t_s[:][0],yerr=avg_t_s[:][1])
plt.ylabel('mean prop distance before <10^12eV or decaying (cm)')
plt.xlabel('energy GeV')
plt.legend(['cont','sto'])
plt.xlim([5,13])
plt.ylim([0,])
plt.show()


