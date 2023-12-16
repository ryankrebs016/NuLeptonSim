import matplotlib.pyplot as plt
import numpy as np

sto=np.loadtxt("out_sto.txt")
cont=np.loadtxt("out_cont.txt")
sto_en=np.loadtxt("sto_en.txt")
cont_en=np.loadtxt("cont_en.txt")

print("avg sto is %f and average cont is %f"%(np.mean(sto)/10**5,np.mean(cont)/10**5))
print("avg sto en is %f and average cont en is %f"%(np.log10(np.mean(sto_en)),np.log10(np.mean(cont_en))))
plt.figure(1)
plt.hist(sto/10**5,histtype='step',bins=np.arange(0,100,1))
plt.hist(cont/10**5,histtype='step',bins=np.arange(0,100,1))
plt.title("Tau flux at 1000 EeV")
plt.legend(['Stochastic','Continuous'])
plt.xlabel("Propagation distance (km)")
plt.ylabel("Count")
#plt.xscale('log')

plt.figure(2)
plt.hist(np.log10(sto_en),histtype='step',bins=np.arange(0,12,.2))
plt.hist(np.log10(cont_en),histtype='step',bins=np.arange(0,12,.2))
plt.title("energy at end")
plt.xlabel("log(GeV)")
plt.legend(['sto','cont'])


#plt.show()
