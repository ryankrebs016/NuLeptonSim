import numpy as np
#angles = 90+np.array([0.1,0.2,0.3,0.5,0.7,0.9,1.2,1.6,2.0,2.4,2.8,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10,12,14,16,18,20,25,30])
#angles = 90+np.array([0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.7,0.9,1.2,1.6,2.0,2.4,2.8,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10,12,14,16,18,20,25,30])
#angles = 90+np.array([0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.05,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10,11,12,13,14,15,16,17,18,19,20,22.5,25,27.5,30])
#angles = 90+np.array([0.1,0.15,0.2,0.25,0.3])
num_parts=1000
energies = np.logspace(3,12,19)
print(energies,np.log10(energies))

#energies = np.sort(np.append(np.logspace(np.log10(2.5e15),np.log10(1.5e16),99),6.3e15))

write_file = open('run_compare_losses',"w")

for j in range(0,len(energies)):
    write_file.write("./loss_comp" +str(" {:.2e} ".format(energies[j])) + str(15)+" "+str(num_parts)+"\n")
    #if j<len(energies)-1: write_file.write("\n")
#write_file.close()

#write_file = open('run_compare_losses_muon', "a")

for j in range(0,len(energies)):
    write_file.write("./loss_comp" +str(" {:.2e} ".format(energies[j])) + str(13)+" "+str(num_parts))
    if j<len(energies)-1: write_file.write("\n")
write_file.close()
