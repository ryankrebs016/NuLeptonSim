import numpy as np
import h5py
import sys
from NuRadioMC.EvtGen.generator import write_events_to_hdf5


if len(sys.argv)<2:
    print("supply input event file")
    exit()

#traj_num,p_thrown,shower_index,vert_x,vert_y,vert_z,tra_x,tra_y,tra_z,nu_prim_flavor,E_nu_prim,part_type,E_part,inel,i_type,had_or_em,sto_index

input_file = sys.argv[1]
data = np.loadtxt(input_file, skiprows=2, unpack=True, delimiter=",")
new_file = input_file.replace(".dat", ".hdf5")
temp = open(input_file,"r")
temp.readline()
temp_str = temp.readline()
temp.close()
temp_str = temp_str.split(",")
rad = float(temp_str[0].split(":")[1])
depth = float(temp_str[1].split(":")[1])
throws = int(temp_str[2].split(":")[1])
#throw_per_traj = int(temp_str[2].split(":")[1])
throws_per_traj = 1

temp_fi = input_file.split("/")[-1]
nu_type_str = temp_fi.split("_")
if len(nu_type_str)==4:
    nu_type = nu_type_str[0]
else:
    nu_type = nu_type_str[0] + "_" + nu_type_str[1]

if nu_type == "nutau":
    nu = 16
elif nu_type == "numu":
    nu = 14
elif nu_type == "nue":
    nu = 12
elif nu_type == "anti_numu":
    nu = -14
elif nu_type == "anti_nutau":
    nu = -16
elif nu_type == "anti_nue":
    nu = -12
else:
    nu=[12, -12, 14, -14, 16, -16]

energy_str = nu_type_str[-2]
if energy_str == "spectrum":
    emin=1e16
    emax=1e21
else:
    emin=float(energy_str)
    emax=float(energy_str)

traj_num = data[0].astype(int)
p_thrown = data[1].astype(int)
show_id = data[2].astype(int)
x = data[3]/1e2
y = data[4]/1e2
z = data[5]/100-6378*1000
tra_x = data[6]
tra_y = data[7]
tra_z = data[8]
prim = data[9].astype(int)
prim_e = data[10] * 1e9
part = data[11].astype(int)
part_e = data[12] * 1e9
inel = data[13]
inter = data[14].astype(int)
shower_type = data[15].astype(int)
sto_i = data[16].astype(int)

nu_i = traj_num+p_thrown*throws_per_traj

r=np.sqrt(tra_x**2+tra_y**2)
azimuth=np.arctan2(tra_y/r,tra_x/r)
zenith=90*np.pi/180-np.arcsin(tra_z)

shower_energy = part_e * inel #eV

attributes = {}
data_sets_fiducial = {}

attributes['NuRadioMC_EvtGen_version'] =-1
attributes['NuRadioMC_EvtGen_version_hash'] = -1
attributes['start_event_id'] = 0
attributes['n_events'] = throws * 1
attributes['flavors'] = nu
attributes['Emin'] = emin
attributes['Emax'] = emax
attributes['thetamin'] = 0
attributes['thetamax'] = np.pi
attributes['phimin'] = 0
attributes['phimax'] = 2*np.pi
attributes['deposited'] = 0

data_sets_fiducial["xx"] = x
data_sets_fiducial["yy"] = y
data_sets_fiducial["zz"] = z
data_sets_fiducial["azimuths"] = azimuth
data_sets_fiducial["zeniths"] = zenith
data_sets_fiducial["event_group_ids"] = nu_i
data_sets_fiducial["flavors"] = part
data_sets_fiducial["energies"] = prim_e
data_sets_fiducial["inelasticity"] = inel
data_sets_fiducial["interaction_type"] = shower_type 
data_sets_fiducial["shower_energies"] = part_e * inel
data_sets_fiducial["weights"] = np.ones(len(shower_type))
temp_show = [""] * len(part)
for i in range(len(temp_show)):
    if shower_type[i] == 0:
        temp_show[i]="had"
    else:
        temp_show[i]="em"
data_sets_fiducial["shower_type"] = temp_show #"had" if shower_type==0 else "em"
#data_sets_fiducial["shower_ids"] = show_id
data_sets_fiducial["shower_ids"] = np.arange(0, len(show_id), 1)


data_sets_fiducial["n_interaction"] = show_id

times = np.zeros(len(sto_i))
for i in range(len(sto_i)):
    if i == 0:
        times[0]=0
    else:
        if nu_i[i] == nu_i[i-1]:
            times[i] = ((z[i]-z[i-1])**2 + (y[i]-y[i-1])**2 + (x[i]-x[i-1])**2)**(1/2) / (2.997925e8*1e-9) + times[i-1]
        else:
            times[i] = 0

data_sets_fiducial["vertex_times"] = times

write_events_to_hdf5(new_file, data_sets_fiducial, attributes, n_events_per_file=throws*throws_per_traj, start_file_id=0)