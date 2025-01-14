import numpy as np
import matplotlib.pyplot as plt
import sys
import math

#vert_x,vert_y,vert_z,
# tra_x,tra_y,tra_z,
# E_nu,inel,part_type,
# i_type,had_or_em,nc_num,
# dc_num,gr_num,traj_num,
# p_thrown,sto_index

dict_data={}
col_map={
    0:"rx",
    1:"ry",
    2:"rz",
    3:"vx",
    4:"vy",
    5:"vz",
    6:"init_flavor",
    7:"E_i",
    8:"E",
    9:"inel",
    10:"part_type",
    11:"int_type",
    12:"shower_type",
    13:"nc_num",
    14:"dc_num",
    15:"gr_num",
    16:"traj_id",
    17:"part_id",
    18:"sto_id"
    }

index=np.arange(0,30,1)

for i in index:
    temp_dict={}

    if len(sys.argv)>2:
        filename=sys.argv[1]
    else:
        filename="testing/events/mixed_events_spectrum_%i.dat"%i

    data=np.loadtxt(filename,skiprows=2,delimiter=",").T

    for ma in range(len(col_map)):
        temp_dict[col_map[ma]]=data[ma]

    if len(dict_data)!=0:
        for key in range(len(col_map)):
            dict_data[col_map[key]]=np.concatenate([dict_data[col_map[key]],temp_dict[col_map[key]]])
    else:
        for key in range(len(col_map)):
            dict_data[col_map[key]]=temp_dict[col_map[key]]

print(dict_data.keys())
print(len(dict_data["rx"]), "events")
dict_data["radius"]=np.sqrt(dict_data["rx"]**2+dict_data["ry"]**2)

dict_data["vxy"]=np.sqrt(dict_data["vx"]**2+dict_data["vy"]**2)

dict_data["fact"]=np.sqrt(dict_data["vx"]**2+dict_data["vy"]**2)
dict_data["azimuth"]=np.arctan2(dict_data["vy"]/dict_data["fact"],dict_data["vx"]/dict_data["fact"])



dict_data["zenith_dir"]=np.arcsin(dict_data["vz"])#np.arctan2(dict_data["vz"],np.sqrt(dict_data["vx"]**2+dict_data["vy"]**2))


plt.figure()
bins=np.arange(11.5,17.5,1)
plt.hist(np.abs(dict_data["part_type"]),bins=bins)
plt.xlabel("part_codes")
plt.xticks([11,12,13,14,15,16])
plt.yscale("log")
plt.savefig("plots/p_types.png")
plt.close()

plt.figure()
bins=np.arange(11.5,17.5,1)
plt.hist(np.abs(dict_data["init_flavor"]),bins=bins)
plt.xlabel("part_codes")
plt.xticks([11,12,13,14,15,16])
plt.savefig("plots/initial_p_types.png")
plt.close()

plt.figure()
bins=np.arange(-.5,7.5,1)
plt.hist(np.abs(dict_data["int_type"]),bins=bins)
plt.xlabel("interaction_types")
plt.xticks([0,1,2,3,4,5,6])
plt.yscale("log")
plt.savefig("plots/interaction_types.png")
plt.close()


plt.figure()
bins=np.arange(-.5,2.6,1)
plt.hist(dict_data["shower_type"],bins=bins)
plt.xticks([0,1,2])
plt.xlabel("shower_types (0=had, 1=EM, 2=HAD+EM)")
plt.savefig("plots/shower_types.png")
plt.close()


plt.figure()
plt.hist(dict_data["radius"]/1e5)
plt.xlabel("vertex radius [km]")
plt.savefig("plots/radius.png")
plt.close()

plt.figure()
plt.hist(dict_data["rx"]/1e5)
plt.xlabel("x [km]")
plt.savefig("plots/rx.png")
plt.close()


plt.figure()
plt.hist(dict_data["ry"]/1e5)
plt.xlabel("y [km]")
plt.savefig("plots/ry.png")
plt.close()

plt.figure()
plt.hist(6378-dict_data["rz"]/1e5)
plt.xlabel("vertex depth [km]")
plt.savefig("plots/depth.png")
plt.close()

plt.figure()
plt.hist(dict_data["vx"])
plt.xlabel("vx")
plt.savefig("plots/vx.png")
plt.close()


plt.figure()
plt.hist(dict_data["vy"])
plt.xlabel("vy")
plt.savefig("plots/vy.png")
plt.close()

plt.figure()
plt.hist(np.log10(dict_data["E"]))
plt.xlabel("particle energy")
plt.savefig("plots/particle_energy.png")
plt.close()


plt.figure()
plt.hist(np.log10(dict_data["E_i"]))
plt.xlabel("primary particle energy")
plt.savefig("plots/initial_primary_energy.png")
plt.close()


plt.figure()
plt.hist(dict_data["zenith_dir"]*180/np.pi)
plt.xlabel("direction elevation angle [deg]")
plt.savefig("plots/elevation_angle.png")
plt.close()


plt.figure()
plt.hist(dict_data["azimuth"]*180/np.pi)
plt.xlabel("Direction Azimuth [deg]")
plt.savefig("plots/azimuth.png")
plt.close()


plt.figure()
plt.hist(np.log10(dict_data["E"]*dict_data["inel"]+.1))
plt.xlabel("log(shower energy [GeV])")
plt.savefig("plots/shower_energy.png")
plt.close()

