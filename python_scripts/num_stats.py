import os
import matplotlib.pyplot as mpl    
    
filename='testing/particles_17.0_99.0.dat'  
print(os.path.exists(filename) )
lc=0    
num_type=[]
for line in open(filename,mode='r',encoding='utf-8'):
    if(lc!=0 and lc!=1 and 'END' not in line):
        num_type.append(int(line.split()[0]))
    lc+=1
muon=0
tau=0
numu=0
nutau=0
for i in range(0,len(num_type)):

    if num_type[i]==1:numu+=1
    if num_type[i]==2:nutau+=1
    if num_type[i]==4:muon+=1
    if num_type[i]==5:tau+=1
print('tau count: %d , muon count: %d , nutau count: %d , numu count: %d '%(tau,muon,nutau,numu))