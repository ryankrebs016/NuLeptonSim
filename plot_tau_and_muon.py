import numpy as np
import matplotlib.pyplot as plt

#energies=np.arange(6,12.5,1)

energies = np.logspace(3,12,19)

energies=np.log10(energies)
for i in range(len(energies)):
  if not (energies[i]%1 ==0):
    energies[i]=energies[i]-0.5+0.4

print(energies)
#print(np.log10(energies))
#for e in np.log10(energies):
#  print(f"{e:0.1f}")

#str(" {:.2e} ".format(energies[j])) 


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
  file_text=f"{e:0.1f}"

  print(f"{e:0.1f}")
  taus_c=np.loadtxt(data_dir+cont_dist_filename+file_text+"_"+str(15)+".txt")*10**-2
  taus_s=np.loadtxt(data_dir+sto_dist_filename+file_text+"_"+str(15)+".txt")*10**-2
  #print(len(taus_c))
  muons_c=np.loadtxt(data_dir+cont_dist_filename+file_text+"_"+str(13)+".txt")*10**-2
  muons_s=np.loadtxt(data_dir+sto_dist_filename+file_text+"_"+str(13)+".txt")*10**-2

  avg_t_c.append(np.mean(taus_c))
  avg_t_s.append(np.mean(taus_s))
  avg_m_c.append(np.mean(muons_c))
  avg_m_s.append(np.mean(muons_s))
  
  #avg_m_c.append([np.mean(muons_c),np.std(muons_c)])
  #avg_m_s.append([np.mean(muons_s),np.std(muons_s)])


plt.figure()
plt.plot(energies,avg_t_c,linestyle='dashed',color='blue')
plt.plot(energies,avg_t_s,color='blue')
plt.plot(energies,avg_m_c,linestyle='dashed',color='red')
plt.plot(energies,avg_m_s,color='red')
plt.ylabel('mean prop distance (m)')
plt.xlabel('energy GeV')
plt.legend(['$\\tau$ cont','$\\tau$ sto','$\\mu$ cont','$\\mu$ sto'])
plt.xlim([2,13])
#plt.yscale('log')


plt.figure()
plt.plot(energies,avg_t_c,linestyle='dashed',color='blue')
plt.plot(energies,avg_t_s,color='blue')
plt.plot(energies,avg_m_c,linestyle='dashed',color='red')
plt.plot(energies,avg_m_s,color='red')
plt.ylabel('mean prop distance (m)')
plt.xlabel('energy GeV')
plt.legend(['$\\tau$ cont','$\\tau$ sto','$\\mu$ cont','$\\mu$ sto'])
plt.xlim([2,13])
plt.yscale('log')


plt.figure()
plt.plot(energies,avg_t_c,linestyle='dashed',color='blue')
plt.plot(energies,avg_t_s,color='blue')
plt.plot(energies,avg_m_c,linestyle='dashed',color='red')
plt.plot(energies,avg_m_s,color='red')
plt.ylabel('mean prop distance (m)')
plt.xlabel('energy GeV')
plt.legend(['$\\tau$ cont','$\\tau$ sto','$\\mu$ cont','$\\mu$ sto'])
plt.xlim([6,13])
#plt.yscale('log')

plt.figure()
plt.plot(energies,avg_t_c,linestyle='dashed',color='blue')
plt.plot(energies,avg_t_s,color='blue')
plt.plot(energies,avg_m_c,linestyle='dashed',color='red')
plt.plot(energies,avg_m_s,color='red')
plt.ylabel('mean prop distance (m)')
plt.xlabel('energy GeV')
plt.legend(['$\\tau$ cont','$\\tau$ sto','$\\mu$ cont','$\\mu$ sto'])
plt.xlim([6,13])
plt.yscale('log')
plt.show()



