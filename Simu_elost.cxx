//############################################################################# 
// Event generator or particle survival code using NuLeptonSim framework
//
// Adapted from code used to obtain probability of emerging tau leptons in 
// calculations of exposure to Earth-skimming neutrinos in Auger and later
// NuTauSim.
//############################################################################# 
// Processes:
// - CC or NC interaction of nu_tau, nu_mu, nu_e in Earth (variable density along chord) 
// - Production of tau and muon leptons with sampling of (1-y) where E_tau=(1-y)*E_nu_tau 
// - Glashow Resonance for anti E neutrinos leading to lepton - neutrino pairs from W boson decay
// - Propagation of leptons (including energy loss)
// - Tau decay sampled from pythia look up tables for resulting neutrino and lepton-neutrino pairs
// - Reinteraction of nu_tau produced in tau decay 
// - Reinteraction of nu_tau produced in nu_tau NC interaction
//----------------------------------------------------------------------------- 
// Several models of neutrino cross-section & lepon energy loss can be chosen
//----------------------------------------------------------------------------- 
// Energies (GeV), unless otherwise specified.
//----------------------------------------------------------------------------- 
// Config:
// - data_dir - Set ouput directory
// - starting_type - Set starting particle type and anti ness *see types below*
// - anti - +1 for normal matter, -1 for anti matter
// - regen - Bool to consider regeneration and decays (0=false,1=true)
// - conversion - Bool to consider simulation of tau decay products
// - bool to use energy distribution, what angles to simulate over
// - Set threshold energy
// - Set detector type (0 - no detector, so particles are earth emerging; 1 - spherical detector, particles saved upon entering volume
//    2 - cylindrical detector)
// - save_neutrinos - Bool to decide if neutrinos neutrinos are saved
// - save_charged - Bool to decide if charged leptons are saved
// - save_events - Bool to decide if events are saved
// - use_sto_inst_of_cont - Bool (True for stochastic, false for continuous)
// - NOT ADDED sto_in_volume - Bool (if true will use stochastic losses inside a prescribed simulation volume to generate detector events) 
//-----------------------------------------------------------------------------
// calling the code is the form ./Event_gen 1E+20 95.0 1E+4 0 0 4.0 0.92 16
// Command line parameters
// ./Event_gen <Initial Neutrino energy> <Emergence Angle> <Number of neutrinos> <Cross section model> <Cross section model> 
//              <Energy loss model> <ice depth> <ice density> <particle type>
// Particle types follow pythia particle codes
// +=Particle, -=Anti Particle. 15=tau, 16=nutau, 13=muons, 14=numus, 12=nue,11=e

//count nu that stop before det, count nu in volume and interact, count nu in volume and leave without
//per flavor per energy

//updated on 7/2/24 by Ryan Krebs

#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include <sys/time.h>
#include <stack>
#include "math.h"
#include <time.h>

#include "Table.hh"
#include "Earth.hh"
#include "Constantes.hh"

#include "sto_losses.h"
#include "cont_losses.h"
#include "Simu_elost.h"


//#define LUT

using namespace std;

//MYEVT_DEF event;
particle_info_t particle_data;
reaction_tables_t reaction_data;
config_t config;
det_t det;



// Initialize Earth class
// The arguments are water thickness and density. 
// They are initialized to bare rock here but it is re-initialized below.
Earth *terra = new Earth(0.0, 2.6); 

//#############################################################
// Main code
//#############################################################
int main(int argc, char **argv)
{
  double time_start = time(NULL); //set start time

  // load config and charged lepton energy loss classes
  load_config(&config); 
  continuous_loss_prop cont;
  stochastic_lepton_prop sto;
  sto.load_tables(config.min_muon_sto_loss, config.min_tau_sto_loss);
  
  
  make_dirs(config.data_dir);
  //load_geo(); not used, things defined in the config file

  int taus_passed_through = 0;
  int decay_num =  0;
  int gr_counter = 0;
  int produced_muons = 0;

  double tau_x = 1;
  double tau_y = 1;
  double tau_z = 1;

  // stuff for external trajectories
  string in_file = "in_files/59068000000_Example_Trajectories.csv";
  int input_num = config.n_traj;
  if(config.detector == 0) input_num = 1;
  double traj_weight = 1/59068000000;
  input_file *input = new input_file[input_num];
  
  // initializes arrays to hold decay products and populates them from pythia file
  for(int i=0; i<100000; i++)
  {
    for(int j=0; j<6; j++)
    {
      reaction_data.tau_type[i][j] = reaction_data.mu_type[i][j] = 0;
      reaction_data.tau_energy[i][j] = reaction_data.mu_energy[i][j] = 0.0;
    }
  }  
  initialize_reaction(reaction_data.tau_type,reaction_data.mu_type,reaction_data.tau_energy,reaction_data.mu_energy);
  
  
  int type_to_save[5] = {-1,-1,-1,-1,-1}; //max of 5 types since electrons are ignored
  
  //add muons and taus to pareticle type to savce
  if(config.save_charged && config.save_emerging)
  {
    if(config.save_mu) type_to_save[0] = 13;
    if(config.save_tau) type_to_save[1] = 15;
  } 
  //add neutrinos to type to save
  if(config.save_neutrinos && config.save_emerging)
  {
    if(config.save_nue) type_to_save[2] = 12;
    if(config.save_numu) type_to_save[3] = 14;
    if(config.save_nutau) type_to_save[4] = 16;
  }
  
  #ifdef DBG
    for(int i=0;i<5;i++) cout<<type_to_save[i]<<" "; //print which particles will be saved
  #endif

  printf("Lepton Propagation code\n");

  
  //-------------------------------------------------
  // Initialisation of ANIS tables for:
  // (a) tau decay - energy of particles produced in tau decay
  // (b) CC and NC interactions of nu_tau - Bjorken y (using CTEQ5)
  // Note: the nu x-section model (Sarkar,etc...)
  //       can be chosen for the propagation of nu_tau through Earth
  //       However CTEQ5 is always used to sample Bjorken y variable.
  //-------------------------------------------------
  
  FinalTable TauData;
  char  taudata[1000];
  char  tfinalccfile[1000];
  char  tfinalncfile[1000];
  char  tfinalccbarfile[1000];
  char  tfinalncbarfile[1000]; 

  FinalTable MuonData;
  char  muondata[1000];
  char  mfinalccfile[1000];
  char  mfinalncfile[1000];
  char  mfinalccbarfile[1000];
  char  mfinalncbarfile[1000];  
  if(argc<=10){
    (void)strcpy(taudata, "tables/tau_decay_tauola.data");
    (void)strcpy(tfinalccfile, "tables/final_cteq5_cc_nu.data");
    (void)strcpy(tfinalncfile, "tables/final_cteq5_nc_nu.data");
    (void)strcpy(tfinalccbarfile, "tables/final_cteq5_cc_nubar.data");
    (void)strcpy(tfinalncbarfile, "tables/final_cteq5_nc_nubar.data");
    //tau above and muon below
    (void)strcpy(muondata, "tables/nu_mu_samples");
    (void)strcpy(mfinalccfile, "tables/final_cteq5_cc_nu.data");
    (void)strcpy(mfinalncfile, "tables/final_cteq5_nc_nu.data");
    (void)strcpy(mfinalccbarfile, "tables/final_cteq5_cc_nubar.data");
    (void)strcpy(mfinalncbarfile, "tables/final_cteq5_nc_nubar.data");
  }

  int InitTau = TauData.InitTable(taudata);
  FinalTable *tCCFinalData = new FinalTable;
  FinalTable *tNCFinalData = new FinalTable;
  FinalTable *tCCBarFinalData = new FinalTable;
  FinalTable *tNCBarFinalData = new FinalTable;
  
  int InitMuon = MuonData.InitTable(muondata);
  FinalTable *mCCFinalData = new FinalTable;
  FinalTable *mNCFinalData = new FinalTable;
  FinalTable *mCCBarFinalData = new FinalTable;
  FinalTable *mNCBarFinalData = new FinalTable;
  tCCFinalData->InitTable(tfinalccfile);
  tNCFinalData->InitTable(tfinalncfile);
  tCCBarFinalData->InitTable(tfinalccbarfile);
  tNCBarFinalData->InitTable(tfinalncbarfile);  

  mCCFinalData->InitTable(mfinalccfile);
  mNCFinalData->InitTable(mfinalncfile);
  mCCBarFinalData->InitTable(mfinalccbarfile);
  mNCBarFinalData->InitTable(mfinalncbarfile);  

  if (argc>1 && atof(argv[1]) == 0) {
    config.energy_distribution = true;
    printf("Will throw uniformly random x neutrinos energy between log10(E_nu/eV) = 15 and  log21(E_nu/eV)\n");
  }
  
  // Initialize Random number generator.
  struct timeval time_struct;
  gettimeofday(&time_struct,NULL);
  srand((time_struct.tv_sec * 1000) + (time_struct.tv_usec / 1000));


  //initialize things and pull in from command line or from config file
  double angle_time_start = time(NULL);
  double angle = config.default_angle;
  if (argc>2) angle = atof(argv[2]);
  if(config.detector == 0) printf("angle %f\n", angle);

  if(config.detector == 0) set_points_from_angle(input, angle);
  if(config.detector == 1)
  {
    if(config.ext_traj) load_input(input, input_num, in_file);
    //else gen_traj();
  }

  // Threshold in which particle will no longer be propagated
  double Elim_eV = config.energy_threshold; // eV
  double Elim = Elim_eV*1.e-9; // GeV
  
  //-------------------------------------------------
  // Declare several variables
  double rndm; // random number throughout the code
  double depth = 0; //center of detector
  double Lmax;
  double dL = 0.;	// Propagation step along chord length. Initialized to zero but it varies for each step.
  double frac = 1e-3; // Fraction of energy lost by tau in each step. This value stays constant
  double dPdes; // Probability of tau decay in dL
  double traversed_grammage; // Track the grammage traversed upon entering the Earth.
  double Energy_GeV; // Particle energy
  double Bjorken_y; // Bjorken y value for interactions CC & NC
  double finalstatecc[2], finalstatenc[2]; // will contain Bjorken (1-y) and y
  float dens; // local density along trajectory
  
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  // Change following parameters depending on application
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  // Number of neutrinos simulated at each zenith angle
  int tot_evt = config.n_throws;
  if(argc > 3)
  {
    tot_evt = (int)atof(argv[3]);
  }
  #ifdef DBG
    printf("Throwing %i neutrinos\n",tot_evt);
  #endif
  /*
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  cout << "======================================" << endl;
  cout << "Number of neutrinos simulated = " << tot_evt << endl;
  cout << "Type of initial neutrinos = "<<config.starting_type<<endl;
  cout << "Energy = " << atof(argv[1]) << " eV" << endl;
  cout << "Threshold energy = " << Elim*1.e9 << " eV" << endl;
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  */

  //for rounding energies and angle in the names

  double starting_energy = config.default_energy;
  if(argc > 1) starting_energy = (double)atof(argv[1]);
  if(config.energy_distribution) starting_energy=0;

  #ifdef DBG
    printf("starting energy 10^%0.2f\n", log10(starting_energy));
  #endif

  string es_temp;
  string e_temp = "";
  if ((int)starting_energy == 0)
  {
    e_temp = "spectrum";
    es_temp = e_temp;
  }
  else
  {
    e_temp=to_string(log10(starting_energy));
    es_temp="";
    for(int i = 0; i<4; i++) es_temp += e_temp[i];
  }

  string ang_temp = to_string(angle);
  string angs_temp = "";
  int count=0;
  if(angle >= 100) count = 6;
  else count = 5;
  for(int i=0; i<count; i++) angs_temp += ang_temp[i];
    
  //cout<<config.starting_type<<","<<atoi(argv[8])<<endl;
  if(argc>9) config.starting_type=atoi(argv[9]);
  printf("overriding starting type as %i\n",config.starting_type);
  //cout<<angs_temp<<endl;
  //-------------------------------------------------
  string tag = "";
  if(argc > 8)
  {
    tag = argv[8];
  }

  // Output file names using input arguments
  string nameEnergies = "";
  if(config.save_emerging) nameEnergies = make_particle_dir(argc, argv, config.data_dir,  es_temp,angs_temp, &config); 

  string nameEvents = "";
  if(config.save_events) nameEvents = make_event_dir(argc, argv, config.data_dir, es_temp, angs_temp, config.starting_type, tag, &config); 

  string name_out_counts = "";
  tag = tag + "_counts";
  if(config.save_final_part_state) name_out_counts = make_event_dir(argc, argv, config.data_dir, es_temp,angs_temp, config.starting_type, tag, &config); 

  //string nameOutLep="";
  //nameOutLep=make_lepton_dir(argc,argv,config.data_dir,es_temp,angs_temp,config.starting_type);
  
  // Open outfiles and place headers
  ofstream out_counts(name_out_counts.c_str());
  out_counts << "nu_id, intial type, ending type, energy_i [GeV], energy_f [GeV], vx, vy, vz, where it stopped {0=b4, 1=in, 2=after}\n";
  
  ofstream outEnergies(nameEnergies.c_str());
  outEnergies << "type, anti, NC, CC, GR, DC, Gen, InitNuNum, InitNeutrinoType, OutEnergy, InitEnergy, Part_Pos.\n";

  ofstream outEvents(nameEvents.c_str());
  outEvents<<"vert_x,vert_y,vert_z,tra_x,tra_y,tra_z,nu_prim_flavor,E_nu_prim,E_part,inel,part_type,i_type,had_or_em,nc_num,dc_num,gr_num,traj_num,p_thrown,sto_index"<< setprecision(9)<<endl;
  
  string out_trajectories_name = "oopsie.txt";
  ofstream out_trajs(out_trajectories_name.c_str());

  #ifdef DBG
    printf("saving emerging particles to %s\n",nameEnergies);
    printf("saving events to %s\n",nameEvents);
  #endif
  outEvents<<"rad [km]: "<<config.ice_det_rad<<", depth [km]: "<<config.ice_det_depth<<", n_nu_throws per file: "<<input_num<<endl;
  

  //ofstream outLep(nameOutLep.c_str());
  //outLep<< setprecision(9)<<"logEi,f,logEd,x,y,z,vx,vy,vz,sto_type,traj,p_num,p_type,sto_num"<<endl;

  //string outGrammage="lepton_test/grammage_"+es_temp+"_"+to_string(config.starting_type)+".dat";

  //ofstream out_gram(outGrammage.c_str());
  //out_gram<<"t_num,n_num,km,grammage,event"<<endl;
  // Get cross-section mode to use
  int CCmode = config.default_cc;
  if(argc>4) CCmode = atoi(argv[4]);
  
  // Get cross-section mode to use
  int ELOSSmode = config.default_eloss;
  if(argc > 5) ELOSSmode = atoi(argv[5]);
  
  // Set water layer properties (bare rock is default)
  terra->depth_new_layer = config.default_layer_thickness;
  terra->dens_new_layer = config.default_layer_density;
  if(argc > 6)
  {
    terra->depth_new_layer = atof(argv[6]);
    terra->dens_new_layer = atof(argv[7]);
    if (terra->dens_new_layer <= 0 || terra->depth_new_layer < 0) {
      cerr << "ERROR: user specified layer most have densitiy >0 and depth >= 0" << endl;
      return -1;
    }
    //cout << "Outer Layer Thickness " << terra->depth_new_layer   << " km" << endl;
    //cout << "Outer Layer Density   " << terra->dens_new_layer    << " g/cm^3" << endl;
  }
  
  // Get zenith angle from input argument
  double refTheta = angle;
  Lmax=2.*R0*cos(PI*(1.-refTheta/180.));
  
  // Average Earth density
  //cout << "Average Earth density " << mean_dens_chord(refTheta) << " g/cm^3" << endl;
  //cout << endl;
  
  cout << "======================================" << endl;
  //cout << "Emerging Lepton Energy file: " << endl;
  //cout << nameEnergies << endl;
  //cout << endl;


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //                                    Initiate run
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //double maxL=sqrt((start_point[0]-end_point[0])*(start_point[0]-end_point[0])+(start_point[1]-end_point[1])*(start_point[1]-end_point[1]));
  //if(Lmax2<=0.) Lmax2=0.;                 // set to zero  if the value is negative.

  if(Lmax<=0.) Lmax=0.;  

  // misc things
  int new_muons = 0; //count of pair of particles made. +1 for muon +1 for neutirno
  int total_CC = 0;
  int total_NC = 0;
  int total_GR = 0;
  double num_tau_decays = 0;
  double num_muon_decays = 0;
  double distance_loop_num = 0;
  int tau_neutrinos_below = 0;
  int muon_neutrinos_below = 0;
  int muons_below = 0;
  int taus_below = 0;
  bool split_type = false;
  int part_count = 0;
  int initial_flavor = 0;
  int event_count = 0;
  int oopsie = 0;  

  if(config.starting_type == 0) split_type = true;

  //start loop over input trajectories
  for(int i=0; i<input_num; i++)
  { 

    if((i+1)%((int)(.1*input_num))==0) printf("ran %i trajectories for %i events\n",i+1,event_count);

    int num_count=0;
    int where_int=0;
    double sum_grammage = 0.;
    double check_grammage=0;
    double maxL;
    double earth_entrance[3];
    double earth_exit[3];
    
    if(!config.ext_traj && config.detector==1) generate_trajectory(earth_entrance, earth_exit, config.ice_det_rad, config.ice_det_depth, config.ang_cutoff);
    //cout<<"gen points "<<earth_entrance[0]<<" "<<earth_entrance[1]<<" "<<earth_entrance[2]<<" "<<earth_exit[0]<<" "<<earth_exit[1]<<" "<<earth_exit[2]<<endl;
    else if(config.ext_traj && config.detector==1)
    {
      earth_entrance[0] = input[i].xi;
      earth_entrance[1] = input[i].yi;
      earth_entrance[2] = input[i].zi;

      earth_exit[0] = input[i].eex;
      earth_exit[1] = input[i].eey;
      earth_exit[2] = input[i].eez;
    }
    else //(config.detector==0)
    {
      earth_entrance[0] = input[0].xi;
      earth_entrance[1] = input[0].yi;
      earth_entrance[2] = input[0].zi;
      earth_exit[0] = input[0].eex;
      earth_exit[1] = input[0].eey;
      earth_exit[2] = input[0].eez;
    }  

    if (config.detector == 1) 
      maxL = sqrt((earth_exit[0]-earth_entrance[0])*(earth_exit[0]-earth_entrance[0])
            +(earth_exit[1]-earth_entrance[1])*(earth_exit[1]-earth_entrance[1])
            +(earth_exit[2]-earth_entrance[2])*(earth_exit[2]-earth_entrance[2]));
    else if(config.detector == 0) maxL = Lmax;
    else maxL = 0;    

    //printf("earth entrance %f %f %f\n",earth_entrance[0],earth_entrance[1],earth_entrance[2]);
    //printf("earth exit %f %f %f\n",earth_exit[0],earth_exit[1],earth_exit[2]);
    //if (config.detector==1)maxL=sqrt((input[i].eex-input[i].xi)*(input[i].eex-input[i].xi)
    //                +(input[i].eey-input[i].yi)*(input[i].eey-input[i].yi)
    //                +(input[i].eez-input[i].zi)*(input[i].eez-input[i].zi));
    //cout<<maxL<<endl;
    if(maxL<10)
    {
      printf("max length, %f cm, error somewhere (check trajectory and detector type) - quitting",maxL);
      exit(-1);
    }

    // get direction unit vectors
    double x_step = (earth_exit[0]-earth_entrance[0])/maxL;
    double y_step = (earth_exit[1]-earth_entrance[1])/maxL;
    double z_step = (earth_exit[2]-earth_entrance[2])/maxL;
    double step_dir[3] = {x_step, y_step, z_step};
    
    double temp_pos[3]={earth_entrance[0],earth_entrance[1],earth_entrance[2]};
    //printf("sim_dir: %f, %f, %f\n",x_step,y_step,z_step);
    
    double* cumulative_grammage = new double[1000000];
    double* grammage_distance = new double[1000000]; // in cm
    double d_grammage=0;
    



    #ifdef LUT
      //printf("using grammage lut\n");
      for (int ii=1; ii<=1000000; ii++)
      {
        double dl = maxL/1000000;
        double x_val = maxL * double(ii) / 1000000.;
        
        temp_pos[0] = temp_pos[0]+dl*x_step;
        temp_pos[1]=temp_pos[1]+dl*y_step;
        temp_pos[2]=temp_pos[2]+dl*z_step;
        sum_grammage +=dl*get_dens_from_coords(temp_pos);
        //sum_grammage +=  dx*earthdens(&x_val, &Lmax);
      }
      d_grammage = sum_grammage/1000000.; // g/cm^2

      cumulative_grammage[0] = 0.;
      grammage_distance[0]   = 0.;
      temp_pos[0]=earth_entrance[0];
      temp_pos[1]=earth_entrance[1];
      temp_pos[2]=earth_entrance[2];
      for (int ii=1; ii<=1000000; ii++)
      {
        double dl = d_grammage/get_dens_from_coords(temp_pos);//chnage to lmax2 for icecube
        double l_val = grammage_distance[ii-1];
        cumulative_grammage[ii] = cumulative_grammage[ii-1] + dl*get_dens_from_coords(temp_pos);//change to lamx2 for icecube
        grammage_distance[ii] = l_val + dl;  
        temp_pos[0] = temp_pos[0]+dl*x_step;
        temp_pos[1] = temp_pos[1]+dl*y_step;
        temp_pos[2] = temp_pos[2]+dl*z_step;
        //printf("*** ii %d %1.5f %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
      }
    #endif

    //----------------------------------
    //cout<<"primary trajectory: "<<i<<endl;
    bool got_event=false;//!got_event && num_count<1000 line below
    int out_leptons=0;


    while((config.run_throws && (num_count<config.n_throws)) || (config.run_number && (out_leptons<config.num_emerging_leptons)))
    //while(out_leptons<500)
    {
      //if(num_count%10000==0) printf("did %i thrown particles\n",num_count);
      //if(out_leptons%100==0) printf("recored %i emerging leptons\n",out_leptons);

      num_count++;

      double pos[3]={earth_entrance[0], earth_entrance[1], earth_entrance[2]};

      // prop specific checks
      bool has_been_saved=false;
      bool save_cond=false;
      bool exit_cond=false;

      // get this particles energy
      Energy_GeV = starting_energy*pow(10,-9);
      if (config.energy_distribution) Energy_GeV =  pow(10,6 + (6 * (double)rand()/RAND_MAX));

      // get this particles type
      int starting_type = 0;
      int part_type = 0;
      int anti = 1;
      anti = config.anti;

      if(config.starting_type==16) 
      {
        part_type=16;
        anti=1;
      }
      else if(config.starting_type==-16)
      {
        part_type=16;
        anti=-1;
      }
      else if(config.starting_type==14)
      {
        part_type=14;
        anti=1;
      }
      else if(config.starting_type==-14)
      {
        part_type=14;
        anti=-1;
      }
      else if(config.starting_type==12)
      {
        part_type=12;
        anti=1;
      }
      else if(config.starting_type==-12)
      {
        part_type=12;
        anti=-1;
      }
      else if(config.starting_type==0)
      {
        part_type=0;
        anti=1;
      }
      else
      {
        printf("use primary neutrino types (12, -12, 14, -14, 16, -16, 0)\n");
        exit(-1);
      }

      //overrule the primary type if split type is true
      if(split_type == true)
      {
        // randomly get flavor
        double rand_type = (double)rand()/RAND_MAX;
        if(rand_type < (1./3.)) part_type = 12;
        else if(rand_type >= (1./3.) && rand_type < (2./3.)) part_type = 14;
        else part_type = 16;
        
        // randomly get if anti particle
        double rand_antiness = (double)rand()/RAND_MAX;
        if(rand_antiness <= (1./2.)) anti = 1;
        else anti = -1;

        #ifdef DBG
          printf("override starting type to use equal distribution of primary flavors and anti-ness\n");
          printf("rand num for type %f\n",rand_type);
          printf("rand num for antiness %f\n",rand_antiness);
          printf("type: %i, antiness: %i\n\n",part_type,anti);
        #endif
      }

      // make internal type
      starting_type = anti*part_type;
  
      if(part_type==0)
      {
        cout<<"Particle type failed to intialize...breaking code\n";
        exit(1);
      }

      initial_flavor = part_type*anti;
      part_count++;

      // various things to track for this particle
      double part_energy = Energy_GeV;
      double part_pos = 0;
      int generation = 0;
      int NC_num = 0;
      int CC_num = 0;
      int GR_num = 0;
      int dc_num = 0;
      
      if(part_pos > maxL || (pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]) > R02)
      {
        // bug catching
        printf("%f,%f,    %f,%f,%f   %f,%f\n",part_pos,maxL,pos[0],pos[1],pos[2],pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2],R02);
        oopsie+=1;
        out_trajs << pos[0] << ","<< pos[1] <<"," <<pos[2] <<","<< x_step << ","<< y_step <<","<< z_step <<"\n";
      }


    
      traversed_grammage = 0.; //initialize traversed grammage for each event
      double traversed_path = 0.;
      
      // Loop over extra particle stacks until empty
      do 
      { 
        // reset some tracking things
        has_been_saved=false;
        bool save_cond=false;
        bool exit_cond=false;
        int start_in_volume=0;
        bool save_lep=false;

        // pop all particle data from the stacks to working variables
        if(!particle_data.part_type.empty())
        {
          part_type = particle_data.part_type.top();
          part_energy = particle_data.part_energy.top();
          part_pos = particle_data.part_pos.top();
          anti = particle_data.anti.top();
          generation = particle_data.generation.top();
          NC_num = particle_data.NC_num.top();
          CC_num = particle_data.CC_num.top();
          dc_num = particle_data.dc_num.top();
          GR_num = particle_data.GR_num.top();
          pos[0] = particle_data.x_pos.top();
          pos[1] = particle_data.y_pos.top();
          pos[2] = particle_data.z_pos.top();
          traversed_grammage = particle_data.traversed_gram.top();

          if(config.save_sto_events) start_in_volume = particle_data.start_in_volume.top();

          particle_data.part_type.pop();
          particle_data.GR_num.pop();
          particle_data.part_energy.pop();
          particle_data.part_pos.pop();
          particle_data.anti.pop();
          particle_data.generation.pop();
          particle_data.NC_num.pop(); 
          particle_data.CC_num.pop();
          particle_data.dc_num.pop();
          particle_data.x_pos.pop();
          particle_data.y_pos.pop();
          particle_data.z_pos.pop();
          particle_data.traversed_gram.pop();

          if(config.save_sto_events) particle_data.start_in_volume.pop();
        }

        double tau_path=0;
        
        // check to make sure electrons or particle below threshold made it into the stack
        if(part_type==11) continue; //ignores particles of electron flavor
        if(part_energy<Elim) continue; //ignore particles below threshold in case they make it through

        bool checked_in=false;
        int sto_num=0;      
        bool broken=false;
        bool left_volume=false;
        bool entered_volume=false;
        int sto_index=0;

        // core propagation loop
        while(part_pos < maxL && (pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]) < R02)
        {
          bool in_vol = in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth);
          
          //create holding arrays for reactions
          int reaction_types[6] = {0,0,0,0,0,0};
          double reaction_energies[6] = {0,0,0,0,0,0};
          int anti_type[6] = {1,1,1,1,1,1};

          // Get the local density for this pa rt of the chord.
          dens = get_dens_from_coords(pos,terra);
          if(dens<0) printf("local dens is < 0\n");

          bool change=false;

          //===========================
          // Particle is a neutrino
          //===========================

          if(part_type == 12 || part_type == 14 || part_type == 16) 
          {
            // Number of interaction lengths propagated in this step is given by an exponentially distributed random number.
            double num_int_lens = -log((double) rand() / (double)(RAND_MAX)); // Randomly sampled number of interaction lengths.
            
            // Convert the number of interaction lengths to grammage.
            double X_int;
            X_int = num_int_lens /(Navo*(dsigGR(part_energy,part_type,anti)+dsigCC(part_energy, CCmode,part_type,anti)+dsigNC(part_energy, CCmode,part_type, anti)));
            
            traversed_grammage += X_int; // Add X_int to the total grammage traversed by the particle
            double nu_step_length = 0; // Initialize the interaction length distance for this step to zero.

            #ifndef LUT      
              // If not using a grammage LUT, a modified earth model is used to analytically calculate step lengths      
              double step_dir[3] = {x_step,y_step,z_step};
              nu_step_length = analytic_distance_from_grammage(pos,step_dir,X_int,terra->depth_new_layer,terra->dens_new_layer);
              part_pos += nu_step_length;
              pos[0] = pos[0]+nu_step_length*x_step;
              pos[1] = pos[1]+nu_step_length*y_step;
              pos[2] = pos[2]+nu_step_length*z_step;
            #endif

            #ifdef LUT
              // The following lines use the grammage_distance lookup table to estimate what position along the trajectory this Xint corresponds to.
              if(traversed_grammage/d_grammage + 1. > 1000000.)
              {
                part_pos = maxL+1;
                traversed_grammage = sum_grammage;
                pos[0] = earth_exit[0];
                pos[1] = earth_exit[1];
                pos[2] = earth_exit[2];  
              }

              // If contained within the trajectory, linearly interpolate its interaction distance.
              if ( floor(traversed_grammage/d_grammage) + 1. < 1000000.) // NOTE: 1000000. is the size of the look-up table.
              {
                double before_step = part_pos;

                // Get the entry in the look-up table corresponding to the traversed grammage
                int ii_grammage = int(traversed_grammage/d_grammage) + 1;
          
                // Linearly interpolate to estimate the distance propagated
                double slope = (grammage_distance[ii_grammage] - grammage_distance[ii_grammage-1])/d_grammage;
                double intercept = grammage_distance[ii_grammage] - slope*cumulative_grammage[ii_grammage];
                //Lint = slope*traversed_grammage + intercept - part_pos; // keep track of this step's interaction length.
                part_pos = slope*traversed_grammage + intercept ;
                double after_step = part_pos;
                pos[0] = pos[0]+(after_step-before_step)*x_step;
                pos[1] = pos[1]+(after_step-before_step)*y_step;
                pos[2] = pos[2]+(after_step-before_step)*z_step;
              }
            #endif
            
            //if(in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth)&& !in_vol)
            //{
            //  entered_volume=true;
              
            //}

            // if the neutrino interaction is still inside Earth, simulate a NC, CC, or GR interaction and check that particle is still above the tracking energy threshold (Elim)
            if(part_pos<maxL)
            {
              double this_dsigNC = dsigNC(part_energy, CCmode, part_type, anti);
              double this_dsigCC = dsigCC(part_energy, CCmode, part_type, anti);
              double this_dsigGR = dsigGR(part_energy, part_type, anti);
              double total_cross = this_dsigNC + this_dsigCC + this_dsigGR;

              double rand_val=((double) rand() / (double)(RAND_MAX));
              if(total_cross==0) printf("total cross section is zero\n");

              double CCratio = this_dsigCC/total_cross;
              double NCratio = this_dsigNC/total_cross;
              double GRratio = this_dsigGR/total_cross;

              bool CChappens = rand_val < CCratio;
              bool NChappens = (CCratio <= rand_val && rand_val < CCratio+NCratio);
              bool GRhappens = (CCratio+NCratio <= rand_val && rand_val < CCratio+NCratio+GRratio);

              if(CChappens)
              {
                //=======================
                // CC interaction occurs (the tracked particle changes from tau neutrino to tau lepton.)
                //=======================
                //total_CC++;
                //cout<<"CC-";
                //CC_num++;
                //cout<<"charged current of neutrinos happened \n";
                // Obtain Bjorken y
                if(part_type == 16)
                {
                  // nutau
                  if(anti == 1) tCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                  if(anti == -1) tCCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                  tau_x = pos[0];
                  tau_y = pos[1];
                  tau_z = pos[2];
                }
                if(part_type==14||part_type==12)
                {
                  //numu or nue
                  if(anti == 1) mCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                  if(anti == -1) mCCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                  if(part_type == 12) broken = true;
                }

                Bjorken_y=finalstatecc[1];
              
                // Set the charged lepton energy from the sampled Bjorken y.
                double initial_energy = part_energy;
                part_energy = (1.-Bjorken_y)*part_energy;
                double shower_energy = initial_energy-part_energy;

                if(shower_energy>Elim && config.save_events && config.save_nu_events && in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
                {
                  // save the event if in volume, part of the events to save, and above threshold
                  got_event = true;
                  event_count++;

                  int shower_code=0;
                  if(part_type == 12)
                  {
                    Bjorken_y = 1;
                    shower_code = 2;
                  }

                  //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
                  outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<< Energy_GeV << ","<< initial_energy<<","<<Bjorken_y<<","
                  <<part_type*anti<<","<<0<<","<<shower_code<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<"\n";

                }

                // Increment the particle counter in the event structure
                //event.npart++;
                //npart++;
                generation++;
                
                if(part_type == 12) part_type = 11; //nue -> e
                if(part_type == 14) part_type = 13; //nu_mu -> mu
                if(part_type == 16) part_type = 15; //nu_tau -> tau

                change=true;
                CC_num++;

              }
              else if (NChappens)
              {
                //=======================
                // NC interaction occurs (the tracked particle remains a tau neutrino with reduced energy.)
                //=======================
                //total_NC++;
                //cout<<"neutral current of neutrinos happened \n";
                // Obtain Bjorken y
                //cout<<"NC-";
                
                if(part_type == 16)
                {
                  if(anti == 1) tNCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                  if(anti == -1) tNCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                
                }
                if(part_type == 14 || part_type == 12)
                {
                  if(anti == 1) mNCFinalData->ThrowFinal(log10(part_energy),finalstatenc);
                  if(anti == -1) mNCBarFinalData->ThrowFinal(log10(part_energy),finalstatenc);
                }
                
                Bjorken_y = finalstatenc[1]; 
                
                // Set the neutrino energy from the sampled Bjorken y.
                double initial_energy = part_energy;
                part_energy = (1.-Bjorken_y)*part_energy;
                double shower_energy = initial_energy-part_energy;

                generation++;

                if(shower_energy>Elim && config.save_events && config.save_nu_events && in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
                {
                  //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
                  outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<<Energy_GeV << ","<<initial_energy<<","<<Bjorken_y<<","
                      <<part_type*anti<<","<<1<<","<<0<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<"\n";
                  event_count++;
                  got_event=true;
                }
                NC_num++;
                change=false;
              }
              else if (GRhappens)
              {
                //total_GR++;
                //GR_num++;
                double initial_energy = part_energy;
                //cout<<"GR-";

                if(part_type!=12 || anti!=-1) 
                {
                  // Check if not a anti nue
                  cout<<CCratio<<" "<<NCratio<<" "<<GRratio<<endl;
                  cout<<"CC happens = "<<CChappens<<" NChappens = "<<NChappens<<" GR happens = "<<GRhappens<<endl;
                  cout<<"GR without a valid particle"<<endl;
                  continue;
                }
                
                double react_rand = ((double) rand() / (double)(RAND_MAX));
                //cout<<react_rand<<endl;
                int temp_channel = -1;
                if(react_rand < 0.676)
                { 
                  temp_channel = 0;
                  //cout<<"W+ decayed to quarks"<<endl;
                  //add in config.save_nu_ev
                  if(initial_energy>Elim && config.save_events && config.save_nu_events && in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
                  {
                    // save hadronic shower
                    //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
                    outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<<Energy_GeV << ","<<initial_energy<<","<<1<<","
                        <<part_type*anti<<","<<2<<","<<temp_channel<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<"\n";
                    
                    event_count++;
                    got_event=true;
                  }

                  GR_num++;
                  broken=true;
                }
                else
                {
                  //cout<<"W+ decayed to leptons"<<endl;
                  double lep_rand = ((double) rand() / (double)(RAND_MAX));
                  double initial_E = part_energy;
                  double lepton_mass = 0;
                  double lepton_type = 0;
                  int lepton_anti = 1;
                  bool has_shower = false;
                  if(lep_rand < (1./3.))
                  {
                    //e made
                    part_type = 12;
                    lepton_mass = me;
                    lepton_type = 11;
                    temp_channel = 1;
                    has_shower = true;
                  }
                  else if(lep_rand < (2./3.) && lep_rand >= (1./3.))
                  {
                    //mu made
                    part_type = 14;
                    lepton_mass = mmuon;
                    lepton_type = 13;
                    temp_channel = 2;
                    has_shower = false;
                  }
                  else if(lep_rand >= (2./3.))
                  {
                    //tau made
                    part_type = 16;
                    lepton_mass = mtau;
                    lepton_type = 15;
                    temp_channel = 3;
                    has_shower = false;
                  }

                  //cout<<"GR happened at "<<pos[0]<<" "<<pos[1]<<" and decayed to "<<part_type<<" and "<<lepton_type<<endl;
                  //cout<<part_type<<","<<lepton_type<<endl;
                  part_energy = ((double) rand() / (double)(RAND_MAX))*part_energy*(mW*mW-lepton_mass*lepton_mass)/(mW*mW);
                  double angle_rand = (double)rand()/(double)RAND_MAX;
                  double lepton_energy = initial_E/mW*(mW*mW+lepton_mass*lepton_mass)/(2*mW)*cbrt(8-8*angle_rand);
                  part_energy = initial_E-lepton_energy;

                  double gr_inel = 0;
                  if(part_type*anti == 12) gr_inel = (initial_E-part_energy)/initial_E;
                  double shower_energy = initial_E*gr_inel;

                  if(has_shower && shower_energy>Elim && config.save_events && config.save_nu_events && in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
                  {
                    //save shower
                    outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<<Energy_GeV << ","<<initial_E<<","<<gr_inel<<","
                        <<part_type*anti<<","<<2<<","<<temp_channel<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<"\n";
                    event_count++;
                    //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
                    got_event=true;
                  }

                  if(lepton_type!=11) 
                  {
                    //throw extra lepton into stack and continue on with the produced nu
                    particle_data.part_type.push(lepton_type);
                    particle_data.part_energy.push(lepton_energy);
                    particle_data.part_pos.push(part_pos);
                    particle_data.anti.push(1);
                    particle_data.generation.push(generation);
                    particle_data.NC_num.push(NC_num);
                    particle_data.CC_num.push(CC_num);
                    particle_data.dc_num.push(dc_num);
                    particle_data.GR_num.push(GR_num);
                    particle_data.x_pos.push(pos[0]);
                    particle_data.y_pos.push(pos[1]);
                    particle_data.z_pos.push(pos[2]);
                    particle_data.traversed_gram.push(traversed_grammage);
                    if(config.save_sto_events) particle_data.start_in_volume.push(in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth));
                  }
                  
                }
                change=true;
                GR_num++;
              }
            } // end of if inside earth block
          }// end of if a neutrino block
          
          //=========================
          // Particle is a tau lepton or muon lepton
          //=========================
          else if(part_type == 11 || part_type == 13 || part_type == 15) //change particle types
          {

            double frac_loss = 0;
            double sampled_decay_length = 0;
            double random_num = 0;

            //if(!in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
            //{
            //cont.set_values(part_energy,dens,part_type,0);
            //dL=cont.get_interaction_length();
            //frac_loss=cont.get_energy_loss()/part_energy;

            //}
            //else 
            //{
            //  sto.set_val(part_energy,part_type,dens);
            //  dL=sto.get_interaction_length();
            //  frac_loss=sto.get_sampled_energy();
            // }
      
            //get a decay length
            random_num = (double)rand()/(double)RAND_MAX;
            sampled_decay_length = decay_length((double)rand()/(double)RAND_MAX,part_energy,part_type); //cm

            // get the interaction info
            // force stochastic if in-ice det and below force distance set in config, or only use stochastic losses
            if((config.detector==1 && (pos[0]*pos[0]+pos[1]*pos[1]+(R0-pos[2])*(R0-pos[2]))<config.sto_force_distance*config.sto_force_distance)||config.use_sto_inst_of_cont)
            {
              sto.set_val(part_energy, part_type, dens);
              dL = sto.get_interaction_length();
              frac_loss = sto.get_sampled_energy();
            }
            else
            {
              // use continuous
              cont.set_values(part_energy, dens, part_type, 0);
              dL = cont.get_interaction_length();
              frac_loss = cont.get_energy_loss()/part_energy;
            }

            //if(tau_path==0)cout<<"tau energy start "<<part_energy<<endl;
            //tau_path+=dL;
            //cout<<tau_path/100000<<"__";
            //cout<<dL<<endl;
            // Check if tau leaves the Earth after dL. If it does then adjust last step
            //if(part_pos+dL > maxL) 
            //{
            //  dL=maxL-part_pos;//change tolmax1 for icecube
            //  //cout<<"tau left before decaying "<<endl;
            //}
            // Calculate the traversed grammage
            //traversed_grammage+=dL*dens;
            

            if(in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth) && !checked_in)
            {
              taus_passed_through++;
              checked_in=true;
            }


            if(sampled_decay_length > dL)
            {
              //==============================
              // The tau lepton does NOT decay
              //=============================
              // cout<<"no decay"<<endl;

              // check to save the deposition
              if(config.save_events && config.save_sto_events && in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth) && ((frac_loss*part_energy)>Elim))
              {
                int temp_show = 0;
                if(sto.sto_type == 0) temp_show=1;
                if(sto.sto_type == 1) temp_show=1;//sto_type==0 brem, sto_type==1 pp, sto_type==2 pn
                if(sto.sto_type == 2) temp_show=0; 
                outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<<Energy_GeV << ","<<part_energy<<","<<frac_loss<<","
                  <<part_type*anti<<","<<sto.sto_type+4<<","<<temp_show<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<sto_num<<"\n";
                //outLep<<log10(part_energy)<<","<<frac_loss<<","<<log10(part_energy*frac_loss)<<","<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<sto.sto_type<<","<<i<<","<<num_count<<","<<part_type<<","<<sto_num<<endl;
                sto_num++;
              }

              // update position and energy
              part_pos = part_pos+dL;
              pos[0] = pos[0]+dL*x_step;
              pos[1] = pos[1]+dL*y_step;
              pos[2] = pos[2]+dL*z_step;
              part_energy=(1.-frac_loss)*part_energy;

              if(part_pos+dL > maxL) 
              {
                dL=maxL-part_pos;
              }

              // Calculate the traversed grammage
              traversed_grammage+=dL*dens;
            }
            else
            {
              //=======================
              // The lepton decays
              //=======================
              // tau or muon does decay
              save_lep=true;
              if(in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
              {
                decay_num++;
              }

              part_pos = part_pos+sampled_decay_length;
              pos[0] = pos[0]+sampled_decay_length*x_step;
              pos[1] = pos[1]+sampled_decay_length*y_step;
              pos[2] = pos[2]+sampled_decay_length*z_step;

              if(part_pos+sampled_decay_length > maxL) 
              {
                sampled_decay_length = maxL-part_pos;//change tolmax1 for icecube

              }
              // Calculate the traversed grammage
              traversed_grammage+=sampled_decay_length*dens;
            
              //dc_num++;
              //cout<<"tau path length before decaying is "<<tau_path<<endl;
              //cout<<"tau energy at decay "<<part_energy<<endl;
              //cout<<"dL is "<<dL<<endl;
              
              // Get the energy of the neutrino produced in the decay
              generation++;

              //if(tag==1) TauData.ThrowFinal(finalstate);
              //if(tag==3) MuonData.ThrowFinal(finalstate);
              //Energy_GeV=finalstate[0]*Energy_GeV;

              double initial_energy = part_energy;
              int reaction_index = (double)rand()/(double)RAND_MAX*100000;

              int initial_particle_type = part_type;
              double dec_inel = 0;
              int initial_particle = part_type;
              int is_had_or_em = -1;
              if(part_type == 15)
              { 
                // is a tau decay
                if(config.conversion)
                {
                  // keep extra particle produced in the decay
                  for(int j=1; j<6; j++)
                  {
                    if(reaction_data.tau_type[reaction_index][j] != 0)
                    { 
                      int anti_ness = 1;
                      if(reaction_data.tau_type[reaction_index][j]<0) anti_ness = -1;
                      anti_type[j] = anti*anti_ness;
                      reaction_types[j] = abs(reaction_data.tau_type[reaction_index][j]);
                      reaction_energies[j] = part_energy*reaction_data.tau_energy[reaction_index][j];
                    }
                    // get what kind of shower will be produced !!!double check
                    if(is_had_or_em == -1 && reaction_types[j] == 0) is_had_or_em = 0;//for HAD
                    if(is_had_or_em == -1 && reaction_types[j] == 11) is_had_or_em = 1;//for EM
                    if(is_had_or_em == -1 && reaction_types[j] == 13) is_had_or_em = 2;//for HAD+EM
                  }
                }
                part_energy *= reaction_data.tau_energy[reaction_index][0];
                part_type = 16;
              }
              if(part_type == 13)
              {
                // is a muon decay
                if(config.conversion)
                {
                for(int j=1; j<6; j++)
                  {
                    if(reaction_data.mu_type[reaction_index][j] != 0) 
                    { 

                      int anti_ness=1;
                      if(reaction_data.mu_type[reaction_index][j]<0) anti_ness = -1;
                      anti_type[j] = anti*anti_ness;
                      reaction_types[j] = abs(reaction_data.mu_type[reaction_index][j]);
                      reaction_energies[j] = part_energy*reaction_data.mu_energy[reaction_index][j];
                    }
                    if(is_had_or_em==-1 && reaction_types[j]==0) is_had_or_em=0;//for HAD
                    if(is_had_or_em==-1 && reaction_types[j]==11) is_had_or_em=1;//for EM
                  }
                }
                part_energy *= reaction_data.mu_energy[reaction_index][0];
                part_type = 14; 
              }

              double frac_energy_dumped = 0;
              for(int frac=0; frac<6; frac++)
              {
                if(reaction_data.tau_type[reaction_index][frac] == 0 || reaction_data.tau_type[reaction_index][frac] == 11)
                frac_energy_dumped += reaction_data.tau_energy[reaction_index][frac];
              }
              //double lost_energy=0;
              //for(int i=1;i<6;i++)
              // {
              //  if(reaction_types[i]!=4) continue;
              //  lost_energy+=reaction_energies[i];
              //}
              //double shower_energy=initial_energy-part_energy-lost_energy;
              //add in config.save_dec
              if((initial_energy*frac_energy_dumped>Elim) && config.save_events&&config.save_dec_events && in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
              {
                // save the decay event
                outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<<Energy_GeV << ","<<initial_energy<<","<<frac_energy_dumped<<","
                    <<initial_particle*anti<<","<<3<<","<<is_had_or_em<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<"\n";
                event_count++;
                //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
                got_event=true;
              }
              dc_num++;
              
              if(config.conversion)
              {
                // add new particle info to simulate the decay products
                for( int j=1;j<6;j++)
                {
                    
                  if((reaction_energies[j]>Elim) && (part_pos<maxL) && (reaction_types[j]!=0))
                  {
                    particle_data.part_type.push(reaction_types[j]);
                    particle_data.part_energy.push(reaction_energies[j]);
                    particle_data.part_pos.push(part_pos);
                    particle_data.anti.push(anti_type[j]);
                    particle_data.generation.push(generation);
                    particle_data.NC_num.push(NC_num);
                    particle_data.CC_num.push(CC_num);
                    particle_data.GR_num.push(GR_num);
                    particle_data.dc_num.push(dc_num);
                    particle_data.x_pos.push(pos[0]);
                    particle_data.y_pos.push(pos[1]);
                    particle_data.z_pos.push(pos[2]);
                    particle_data.traversed_gram.push(traversed_grammage);
                    if(config.save_sto_events) particle_data.start_in_volume.push(in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth));
                  }
                }
              }
              // if not regenerating to a neutrino, break the loop
              if(!config.regen) broken=true;
            }
          } // end if particle is a charged lepton

          // double check if it's below threshold
          if(part_energy<Elim) broken=true;

          // double check if it's still traveling towards the in-ice det
          double going_still = still_going_towards_det(pos, step_dir, config.ice_det_rad,config.ice_det_depth);
          if(going_still==2) 
          {
            if (config.save_final_part_state) out_counts << i <<","<< starting_type << "," << part_type*anti<<","<< Energy_GeV <<","<< part_energy<<","<<x_step<<","<<y_step<<","<<z_step<<","<<going_still<<"\n";
            if (config.detector==1) break;
          }
          if((part_type!=12 && part_type!=13 && part_type!=14 && part_type!=15 && part_type!=16) || part_energy<Elim || broken || (pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])>R02)
          {
            if(config.save_final_part_state) out_counts << i <<"," << starting_type << "," << part_type*anti<<","<< Energy_GeV <<","<< part_energy<<","<<x_step<<","<<y_step<<","<<z_step<<","<<going_still<<"\n";
            break;
          }
        } // ends core propagation loop
        
        //=================================================
        // Write energy of emerging tau to output text file
        //=================================================
        for(int j=0; j<5;j++)
        {
          if(config.save_emerging && (part_type==type_to_save[j]) && (part_energy>Elim) && broken==false && has_been_saved==false)
          {
            out_leptons++;
            outEnergies <<part_type<<" "<<anti<<" "<<NC_num << " " << CC_num << " " << GR_num<<" "<<
            dc_num << " " <<generation <<" "<<num_count<< " " << log10(part_energy)+9 << " " << log10(Energy_GeV)+9<<" "<<part_pos<<" "<<tau_x<<" "<<tau_y<<" "<<tau_z<<"\n";
            has_been_saved=true;
            
            if(out_leptons==config.num_emerging_leptons) outEnergies<<num_count<<" initial neutrinos of type "<<config.starting_type<<" at energy "<<Energy_GeV*pow(10,9);
          }
        }
      } while(!particle_data.part_type.empty()); // end of stack loop
    } // end of num particles loop

    delete cumulative_grammage;
    delete grammage_distance;
  } // end of trajectories loop
  
  printf("num bad start points %i\n",oopsie);
  printf("taus decayed: %i, taus passing through (could have decayed): %i",decay_num,taus_passed_through);
  //outEnergies << "END" << endl; // write END in the last line of the text output file. 
  outEnergies.flush();
  //outEvents <<"END"<<endl;
  outEvents.flush();
  out_counts.flush();
  //out_gram<<"END"<<endl;
  //out_gram.flush();
  out_trajs.flush();
  //outLep<<"END"<<endl;
  //outLep.flush();
  double angle_time_end = time(NULL);
  cout<<"total CC "<<total_CC<<". total NC "<<total_NC<<". total GR "<<total_GR<<endl;
  //muon_output<<angle<<" "<<produced_muons<<endl;
  //output<<angle<<" "<<tot_evt<<" " <<angle_time_end-angle_time_start<<" "<<endl;

  /*
  cout<<"total CC "<<total_CC<<". total NC "<<total_NC<<endl;
  cout<<"num of CC "<<total_CC<<endl;
  cout<<"num of NC "<<total_NC<<endl;
  cout<<"num tau decays "<<num_tau_decays<<endl;
  cout<<"num muon decays "<<num_muon_decays<<endl;
  cout<<"num distance loop per part "<<distance_loop_num/100000<<endl;
  cout<<"taus below "<<taus_below<<endl;
  cout<<"tau neutrinos below "<<tau_neutrinos_below<<endl;  
  cout<<"muons below "<<muons_below<<endl;
  cout<<"muon neutrinos below "<<muon_neutrinos_below<<endl;
  */
 
  double time_elapsed = time(NULL)-time_start;
  printf("Elapsed time %f with %i events\n",time_elapsed,event_count);
  
  return 0;
  
}    
