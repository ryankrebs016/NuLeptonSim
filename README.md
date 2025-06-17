# README #

Documentation
##################

The code serves 3 main functionalities that will be outline in the next sections. Run specific settings include cross sectional models, charged lepton energy loss methods, primary energies, parimary flavors, threshold energies, trajectories through earth, outer Earth layer properties. Energies can be set in the command line, config, or from a random log uniform distribution. Cross sectional models are taken on command line or config file. Charged lepton energy loss can be continuous, where the particle loses a set percent of energy and the step is adjusted from that, or can be stochastic, where interaction lengths are randomly sampled and energy losses sampled from Bremsstrahlung, pair production, or photoelectric processes. Some physics options for taus (and techinically muons) exists to allow regeneration (tau->nutau) and decay products (tau -> nutau + stuff) in the simulations. Lastly, the outer layer can be specified with a density and depth, and replaces the outer layers of the PREM Earth density model.

All options that appear in the command line should override what they are set to in the config file.


Ways to use NuLeptonSim
-------------------------
1. The first being for high altitude detectors that need exit probabilities and energy distributions of charged secondaries emerging from Earth. In the config file, set detector=0, save_emerging=1, then save_charged=1 to enable saving charged leptons and further specifying which flavor and/or save_neutrinos=1 and the flavors to save emerging neutrinos. Set save_final_particle_state=0. Directions can either be set with the exit angle in the command line, measured from downgoing trajectories (95 would point slightly above the horizon, or in the config file. The time that NLS runs is determined by how many neutrinos it runs. If run_throws=1, it will throw n_throws*n_trajectories total neutrinos at one angle, so n_trajectories should be set to 1. Alternatively, if run_throws=0 it will run until numm_emerging_leptons is reached (this is exceedingly slow for steep exit angles where exit probabilities are very low). The idea here is to run NuLeptonSim over lists of primary flavors and energies to build lookup tables.

2. The second is for in-ice detectors would detect particle showers initiaited by primaries or secondaries interacting in the ice. This method of simultation, instead of looking at emerging particles from a specific directions, will instead want to throw neutrinos into the detector volume specific and save all events that occur inside of it. In the config file, set detector=1 and then specify the volume of a cylindrical detector through it radius and height. The top of the cylinder is centered at (0,0,Re) and exists downward into the ice. Here, save_emerging=0, save_charged=0, save_neutrinos=0, save_final_particle_state=0. Instead use save_events=1, with the subtype of events neutrino events (CC, NC, GR), save dec(ay) events for taus and muons, save sto(chastic) events for taus and muons too. Stochastic events, if set as they do improve detector performance, with dramatically exceed the number of primary neutrinos events if the threshold is below 10PeV - this will slow down the detector side of the simulations so set a reasonable energy threshold. The trajectories neutrinos take are determined by first generating a random point inside the volume, randomly getting a directions, and then backtracking to the surface of Earth as the starting position. The last thing to mention, is that the number of trajectories generated will match n_trajectories in the config and will throw n_throws along each trajectory. So n_throws should be set to 1.

3. This last mode is linked to the in-ice mode. Here instead of saving emerging particles from earth or saving the events in a volume, this mode will save the survival probabilites to the volume. If using events generated from NuLeptonSim, getting an effective volume is somewhat roundabout. First you generate the events in 2, then you need to have dedicated 3. runs to get the probability a particle makes it to the volume and then the probability it interacts with some energy is encoded in the event numbers. One can get an effective area instead, though this requires external trajectories with random start point and end points - NuLeptonSim does not do this as it requires it's own Monte Carlo to generate these trajectores. To run it in this mode, set save_events=0, save_emerging=0, and instead have save_final_particle_state=1. Settings should match what was done for 2 except num_throws can be much bigger as only a few files are needed to be made.



Physics
---------------------
This project is built upon the base code of NuTauSim (https://github.com/harmscho/NuTauSim).

NuTauSim:
The underlying principle of NuTauSim, as well as some studies, are documented in https://arxiv.org/abs/1707.00334  and https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.023021 . NuTauSim can be cited with https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.023021 . A bug fix ( https://github.com/harmscho/NuTauSim/commit/ddacd1aeebab7cdfa6c155464b79d549e8b9b02e ) lead to an erratum which corresponds to the current version of the this code. The erratum you can find here https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.069902 .

NuLeptonSim expanded on NuTauSim by simulating the muons flavors, tau/muon decay products, stochastic energy losses, and tailoring the code to be used with in-ice detectors.
Some work using this code can be found in arXiv:2311.03646, outlining the muon behavior and the stochastic energy losses, and arXiv:2308.07401, which passes NuLeptonSim events to a detector simulation of ARA to estimate the effects of secondaries.

How to compile and run
-----------------------
Step 0:
If using stochastic losses, first enter the stochastic_tables directory and run "setup_tables.sh" to unpack all the tables.

Step 1:
Fromt the main NuLeptonSim directory run "make" to generate executables (binaries)

Step 2:
Edit the config file to have the correct directory to save all the data files and set other settings for how you want to run NuLeptonSim

Step 3:
Run NuLeptonSim from the command line, see below. The command line arg order may change format, so tracking down what arg number is which thing in the source file may be necessary.

As of now, it is called like this:
Simu_elost [Energy eV] [Angle] [Number of Events] [Nu Cross Section Model] [Tau Energy Loss Model] [Water Layer Thickness] [Water Layer Density] [Output Tag] [Primary particle types] 

Example: 
./Simu_elost 1e20 91.0 1e2 0 0 4.0 0.92 test_tmp_or_index 16

This injects tau neutrinos (16, nue 12, numu 14) with energy 10^20 eV at an angle of 91 degrees (exit angle 89 degrees in zenith, emergence angle 1 degree above horizon). 100 (1e2) tau neutrinos are simulated. The cross-section model is the middle (a.k.a. standard) curve (0) and the tau energy loss rate model is ALLM (0) (unless overidden with stochastic loss). The water layer thickness is 4.0 km with density 0.92 g/cm^3 (ice). The output files have tag "test_tmp_or_index" and will be directing into the data_dir in the config.

One can give the energy a value of 0, in this case tau-energies are generated (randomly) uniformly in log-space in the range E=1e15 eV and E=1e21.

Step 4:.
Generate relevant lookup tables. The look-up tables are compiled with make_LUT.py and placed in the LUTs subdirectory in the larger data_dir.
NOTE: the user will have to edit the script to specify their relevant outdir and exe directories. Keep an eye out for any directories written in by hand.
