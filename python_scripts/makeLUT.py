import matplotlib.pyplot as plt
import numpy as np
import os

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
    nu_type=[]
    n_nu=0
    grab_last=False

    for line in open(filename,mode='r',encoding='utf-8'):
        if "Threw" in line:
            n_nu = int(line.split()[1])
        elif lc==1 and "until" not in line:
            n_nu = int(line.split()[1])
        elif(lc!=0 and lc!=1):
            try:
                num_type.append(int(line.split()[0]))
                anti_type.append(int(line.split()[1]))
                num_NC.append(int(line.split()[2]))
                num_CC.append(int(line.split()[3]))
                num_GR.append(int(line.split()[4]))
                num_decays.append(int(line.split()[5]))
                num_gen.append(int(line.split()[6]))
                nu_type.append(int(line.split()[7]))

                #num_particles.append(int(line.split()[7]))
                energy.append(float(line.split()[9]))
                start_energy.append(float(line.split()[10]))
                end_pos.append(float(line.split()[11]))
                #pos=10
            except:
                print("failed", line)
                exit()

        lc+=1
    return np.array(num_type),np.array(anti_type),np.array(num_CC), np.array(num_NC), np.array(num_decays), np.array(num_particles), np.array(energy), n_nu

def map_flavor(num):
    if num==16:
        return "nutau"
    elif num==-16:
        return "anti_nutau"
    elif num==14:
        return "numu"
    elif num==-14:
        return "anti_numu"
    elif num==12:
        return "nue"
    elif num==-12:
        return "anti_nue"
    else:
        raise(ValueError,"incorrect particle code")

data_dir = 'output/'
particle_dir=data_dir+'particles/'
paths=os.listdir(particle_dir)

outdir = data_dir+'LUT/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

angles = np.loadtxt("angles.txt", ndmin=1)
flavors = np.loadtxt("flavors.txt", ndmin=1)
energies = np.loadtxt("energies.txt", ndmin=1)

for flavor in flavors:
    for energy in energies:
        for angle in angles:
            fnm = data_dir+f"particles/{map_flavor(flavor)}_leptons_{np.log10(energy):.1f}_{angle:.2f}.dat"
            the_type,anti_t,num_CC, num_NC, num_decays, num_particles, part_energy, thrown_nu = read_emerging(fnm)

            print(f'P_surv = {len(part_energy)/thrown_nu} for {fnm}' )
            print(part_energy)
            np.savez(f'{outdir}LUT_{map_flavor(flavor)}_leptons_{np.log10(energy):.1f}_eV_{angle:.2f}_deg.npz', n_thrown=thrown_nu,
                     flavor_code=the_type*anti_t, part_energies = part_energy, num_CC=num_CC,num_NC=num_NC,num_decays=num_decays)
