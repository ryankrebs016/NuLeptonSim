import numpy as np
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

import time

#Setting material properties
#Air
#Z, A = 7.42, 15

#Ice
Z, A = 11, 22
rho = 0.92 #g/cm^3

#Loading cdf data table
cdf_data_muon = np.load("Muon_Radiative_Cross_Section_CDFs.npy", allow_pickle = True)
cdf_data_tau = np.load("Tau_Radiative_Cross_Section_CDFs.npy", allow_pickle = True)

MeV_energies_muon, cross_sections_muon, cdf_values_muon, cdf_xs_muon = cdf_data_muon[0], cdf_data_muon[1], cdf_data_muon[2], cdf_data_muon[3]
cross_section_brem_muon, cross_section_pp_muon, cross_section_pn_muon = cross_sections_muon[:,0], cross_sections_muon[:,1], cross_sections_muon[:,2]
cross_section_total_muon = cross_section_brem_muon+cross_section_pp_muon+cross_section_pn_muon
print(np.shape(MeV_energies_muon),np.shape(cross_sections_muon), np.shape(cdf_values_muon),np.shape(cdf_xs_muon))
'''
np.savetxt('stochastic_tables/MeV_energies.txt',MeV_energies_muon)

np.savetxt('stochastic_tables/cs_brem_muon.txt',cross_section_brem_muon)
np.savetxt('stochastic_tables/cs_pp_muon.txt',cross_section_pp_muon)
np.savetxt('stochastic_tables/cs_pn_muon.txt',cross_section_pn_muon)

np.savetxt('stochastic_tables/cdf_values_muon.txt',cdf_values_muon)
np.savetxt('stochastic_tables/cdf_xs_brem_muon.txt',cdf_xs_muon[:,0,:],delimiter=' ')
np.savetxt('stochastic_tables/cdf_xs_pp_muon.txt',cdf_xs_muon[:,1,:],delimiter=' ')
np.savetxt('stochastic_tables/cdf_xs_pn_muon.txt',cdf_xs_muon[:,2,:],delimiter=' ')

'''

MeV_energies_tau, cross_sections_tau, cdf_values_tau, cdf_xs_tau = cdf_data_tau[0], cdf_data_tau[1], cdf_data_tau[2], cdf_data_tau[3]
cross_section_brem_tau, cross_section_pp_tau, cross_section_pn_tau = cross_sections_tau[:,0], cross_sections_tau[:,1], cross_sections_tau[:,2]
cross_section_total_tau = cross_section_brem_tau+cross_section_pp_tau+cross_section_pn_tau
'''
np.savetxt('stochastic_tables/cs_brem_tau.txt',cross_section_brem_tau)
np.savetxt('stochastic_tables/cs_pp_tau.txt',cross_section_pp_tau)
np.savetxt('stochastic_tables/cs_pn_tau.txt',cross_section_pn_tau)

np.savetxt('stochastic_tables/cdf_values_tau.txt',cdf_values_tau)
np.savetxt('stochastic_tables/cdf_xs_brem_tau.txt',cdf_xs_tau[:,0,:],delimiter=' ')
np.savetxt('stochastic_tables/cdf_xs_pp_tau.txt',cdf_xs_tau[:,1,:],delimiter=' ')
np.savetxt('stochastic_tables/cdf_xs_pn_tau.txt',cdf_xs_tau[:,2,:],delimiter=' ')
'''
print(MeV_energies_muon,MeV_energies_tau)

# Interpolating cross sections
cs_brem_muon, cs_pp_muon, cs_pn_muon = interpolate.interp1d(MeV_energies_muon, cross_section_brem_muon), interpolate.interp1d(MeV_energies_muon, cross_section_pp_muon), interpolate.interp1d(MeV_energies_muon, cross_section_pn_muon)
cs_brem_tau, cs_pp_tau, cs_pn_tau = interpolate.interp1d(MeV_energies_tau, cross_section_brem_tau), interpolate.interp1d(MeV_energies_tau, cross_section_pp_tau), interpolate.interp1d(MeV_energies_tau, cross_section_pn_tau)

interaction_length_muon = interpolate.interp1d(MeV_energies_muon, A*1.66e-24/cross_section_total_muon)
interaction_length_tau = interpolate.interp1d(MeV_energies_tau, A*1.66e-24/cross_section_total_tau)

# en = np.linspace(5,15,100)
# fig = plt.figure()
# plt.plot(en,interaction_length_muon(en), label = 'muon')
# plt.plot(en,interaction_length_tau(en), label = 'tau')
# plt.legend()
# plt.show()

#Creating regular grid interpolator for energy sampling
interpolator_function_brem_muon = RegularGridInterpolator((MeV_energies_muon, cdf_values_muon), cdf_xs_muon[:,0,:])
interpolator_function_pp_muon = RegularGridInterpolator((MeV_energies_muon, cdf_values_muon), cdf_xs_muon[:,1,:])
interpolator_function_pn_muon = RegularGridInterpolator((MeV_energies_muon, cdf_values_muon), cdf_xs_muon[:,2,:])

interpolator_function_brem_tau = RegularGridInterpolator((MeV_energies_tau, cdf_values_tau), cdf_xs_tau[:,0,:])
interpolator_function_pp_tau = RegularGridInterpolator((MeV_energies_tau, cdf_values_tau), cdf_xs_tau[:,1,:])
interpolator_function_pn_tau = RegularGridInterpolator((MeV_energies_tau, cdf_values_tau), cdf_xs_tau[:,2,:])


# energies = np.linspace(5,15,100)
# fig = plt.figure()
# plt.plot(energies, interaction_length_muon(energies))
# plt.plot(energies, interaction_length_tau(energies))
# plt.yscale('log')
# plt.show()
# samples = np.random.rand(1000000)
# energies = 5*np.ones(len(samples))
# 
# desired_points = np.array([energies.flatten(), samples.flatten()]).T
# 
# tic = time.time()
# blah1 = interpolator_function_brem_muon(desired_points)
# blah2 = interpolator_function_pp_muon(desired_points)
# blah3 = interpolator_function_pn_muon(desired_points)
# # 
# blah4 = interpolator_function_brem_tau(desired_points)
# blah5 = interpolator_function_pp_tau(desired_points)
# blah6 = interpolator_function_pn_tau(desired_points)
# # print(time.time()-tic)
# # 
# # fig = plt.figure()
# # plt.hist(blah1, bins = np.logspace(np.log10(min(blah1)), np.log10(max(blah1)), 100), density = True, histtype = 'step', color = 'green')
# # plt.hist(blah2, bins = np.logspace(np.log10(min(blah1)), np.log10(max(blah1)), 100), density = True, histtype = 'step', color = 'orange')
# # plt.hist(blah3, bins = np.logspace(np.log10(min(blah1)), np.log10(max(blah1)), 100), density = True, histtype = 'step', color = 'blue')
# # 
# plt.hist(blah4, bins = np.logspace(np.log10(min(blah1)), np.log10(max(blah1)), 100), histtype = 'step', color = 'green', label = 'Brem')
# plt.hist(blah5, bins = np.logspace(np.log10(min(blah1)), np.log10(max(blah1)), 100), histtype = 'step', color = 'orange', label = 'PP')
# plt.hist(blah6, bins = np.logspace(np.log10(min(blah1)), np.log10(max(blah1)), 100), histtype = 'step', color = 'blue', label = 'PN')
# 
# plt.xlabel(r'$y = E/E_{\tau}$')
# plt.ylabel(r'N')
# plt.title(r'$E_{\tau} = 10^{11}$'+' eV')
# plt.xscale('log')
# plt.yscale('log')
# plt.legend()
# plt.grid(which = 'both')
# plt.show()

#Loading charged leptons from NuLeptonSim runs
input_file = 'lepton_test/18.0_16.dat'
Ryan_leptons = np.loadtxt(input_file, delimiter=',', skiprows=1)
output_file = 'lepton_test/18.0_16_leptons.csv'

energies_initial = np.log10(1e3*Ryan_leptons[:,6])
energies = np.log10(1e3*Ryan_leptons[:,6])
lepton_type = Ryan_leptons[:,8]
trajectory = Ryan_leptons[:,10]
throw = Ryan_leptons[:,11]
start_pos = np.array([Ryan_leptons[:,0], Ryan_leptons[:,1], Ryan_leptons[:,2]])
end_pos = np.array([Ryan_leptons[:,3], Ryan_leptons[:,4], Ryan_leptons[:,5]])
'''

energies_initial = np.log10(1e3*Ryan_leptons[1,6])
energies = np.log10(1e3*Ryan_leptons[1,6])
lepton_type = Ryan_leptons[1,8]
trajectory = Ryan_leptons[1,10]
throw = Ryan_leptons[1,11]
start_pos = np.array([Ryan_leptons[1,0], Ryan_leptons[1,1], Ryan_leptons[1,2]])
end_pos = np.array([Ryan_leptons[1,3], Ryan_leptons[1,4], Ryan_leptons[1,5]])
'''
dist1 = np.linalg.norm(end_pos-start_pos, axis = 0)
dir_vec = (end_pos-start_pos)/dist1
dist2 = 1e20*np.ones(dist1.shape)

events = np.empty((0,13), dtype = 'object')
threshold = 1e11 #MeV
loops=0
print(np.shape(energies))
while np.any(np.logical_and(dist2>dist1, energies>np.log10(threshold))):

    loops=loops+1
    idx = np.logical_and(dist2>dist1, energies>np.log10(threshold))
    print('-----------------------------')
    print('idx is ',np.shape(idx),idx)

    
    #print(np.shape(energies),np.shape(idx))
    #print(sum(idx))
    dist2 = np.copy(dist1)
    num_leptons = sum(idx)
    num_int_lengths = np.random.exponential(size = num_leptons)
    print('\nnum len is ',np.shape(num_int_lengths),print(num_int_lengths))
    int_dist = np.empty(shape = num_leptons)
    
    is_muon = lepton_type[idx] == 13
    is_tau = lepton_type[idx] == 15
    
    int_dist[is_muon] = num_int_lengths[is_muon]*interaction_length_muon(energies[idx][is_muon])/rho
    int_dist[is_tau] = num_int_lengths[is_tau]*interaction_length_tau(energies[idx][is_tau])/rho

    start_pos[:, idx] = start_pos[:, idx]+int_dist*dir_vec[:, idx]
    dist1 = np.linalg.norm(end_pos-start_pos, axis = 0)
    
    prob = np.random.rand(num_leptons)
    
    muon_energies, tau_energies = energies[idx][is_muon], energies[idx][is_tau]
    
    brem_muon, pp_muon, pn_muon = cs_brem_muon(muon_energies), cs_pp_muon(muon_energies), cs_pn_muon(muon_energies)
    brem_tau, pp_tau, pn_tau = cs_brem_tau(tau_energies), cs_pp_tau(tau_energies), cs_pn_tau(tau_energies)
    print(np.shape(brem_muon),'----')
    exit()
    total_muon, total_tau = brem_muon+pp_muon+pn_muon, brem_tau+pp_tau+pn_tau

    int_type = np.empty(prob.shape, dtype="<U10")
    int_type[is_muon] = np.where(prob[is_muon]<brem_muon/total_muon, 'brem', np.where(prob[is_muon]<(brem_muon+pp_muon)/total_muon, 'pp', 'pn'))
    int_type[is_tau] = np.where(prob[is_tau]<brem_tau/total_tau, 'brem', np.where(prob[is_tau]<(brem_tau+pp_tau)/total_tau, 'pp', 'pn'))
    print('\n int type is ',np.shape(int_type),int_type)
    sampled_energy = np.empty(prob.shape)
    sampled_energy[np.logical_and(is_muon, int_type == 'brem')] = interpolator_function_brem_muon(np.array([energies[idx][np.logical_and(is_muon, int_type == 'brem')], np.random.rand(len(energies[idx][np.logical_and(is_muon, int_type == 'brem')]))]).T)
    sampled_energy[np.logical_and(is_muon, int_type == 'pp')] = interpolator_function_pp_muon(np.array([energies[idx][np.logical_and(is_muon, int_type == 'pp')], np.random.rand(len(energies[idx][np.logical_and(is_muon, int_type == 'pp')]))]).T)
    sampled_energy[np.logical_and(is_muon, int_type == 'pn')] = interpolator_function_pn_muon(np.array([energies[idx][np.logical_and(is_muon, int_type == 'pn')], np.random.rand(len(energies[idx][np.logical_and(is_muon, int_type == 'pn')]))]).T)

    sampled_energy[np.logical_and(is_tau, int_type == 'brem')] = interpolator_function_brem_tau(np.array([energies[idx][np.logical_and(is_tau, int_type == 'brem')], np.random.rand(len(energies[idx][np.logical_and(is_tau, int_type == 'brem')]))]).T)
    sampled_energy[np.logical_and(is_tau, int_type == 'pp')] = interpolator_function_pp_tau(np.array([energies[idx][np.logical_and(is_tau, int_type == 'pp')], np.random.rand(len(energies[idx][np.logical_and(is_tau, int_type == 'pp')]))]).T)
    sampled_energy[np.logical_and(is_tau, int_type == 'pn')] = interpolator_function_pn_tau(np.array([energies[idx][np.logical_and(is_tau, int_type == 'pn')], np.random.rand(len(energies[idx][np.logical_and(is_tau, int_type == 'pn')]))]).T)
    print('\n int type is ',np.shape(sampled_energy),sampled_energy)
    dep_energy = sampled_energy*10**energies[idx]
    energies[idx] = np.log10(10**energies[idx]-dep_energy)

    bright_events = np.logical_and(dep_energy>threshold, dist2[idx]>dist1[idx]) 

    pos = start_pos[:, idx][:, bright_events]
    dir = dir_vec[:, idx][:, bright_events]
    if loops >10: exit()
    if np.sum(bright_events)>0:
        #print(trajectory[idx][bright_events], throw[idx][bright_events])
        #a_size=np.size(lepton_type[idx][bright_events])
        #print(a_size)
        #print(np.arange(1,a_size,1))
        events = np.append(events, np.array([energies_initial[idx][bright_events], sampled_energy[bright_events], np.log10(dep_energy[bright_events]), pos[0], pos[1], pos[2], dir[0], dir[1], dir[2], int_type[bright_events], trajectory[idx][bright_events], throw[idx][bright_events],lepton_type[idx][bright_events]], dtype = 'object').T, axis=0)
arrform=' '.join(['%f']*9 + ['%s']+['%d']*2+['%i'])

np.savetxt(output_file, events, fmt=arrform, delimiter=",")

# delta = (10**energies_initial-10**energies)/(10**energies_initial)
# 
# fig = plt.figure()
# plt.hist(delta, bins = np.logspace(np.log10(min(delta)),np.log10(max(delta)), 100))
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$E_{f}/E_{i}$')
# plt.title('Energy Losses in ARA from 100EeV primary '+r'$\mu$'+' neutrinos')
# #plt.title('Energy Losses in ARA from 10EeV primary '+r'$\tau$'+' neutrinos')
# 
# 
# fig = plt.figure()
# data = events[:, 0].astype('float')
# plt.hist(data, bins = np.logspace(np.log10(min(data)), np.log10(max(data)), 30), density = True)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Deposited Energy (MeV)')
# plt.ylabel(r'$dn/dE$')
# plt.title('Depositions in ARA from 100EeV primary '+r'$\mu$'+' neutrinos')
# #plt.title('Depositions in ARA from 10EeV primary '+r'$\tau$'+' neutrinos')
# plt.show()


