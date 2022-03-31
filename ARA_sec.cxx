//############################################################################# 
// Propagation of neutrinos in Earth - flux of emerging leptons
//
// Adapted from code used to obtain probability of emerging tau leptons in 
// calculations of exposure to Earth-skimming neutrinos in Auger
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
// - Set ouput directory !!! need to first create that directory before running and the /particles subdirectory !!!
// - Set starting particle type and anti ness (1-E,2-M,3-T,4-nu_E,5-nu_M,6-nu_T) (-1=anti particle, 1=normal)
// - Bool to consider regeneration and decays (0=false,1=true)
// - bool to use energy distribution, what angles to simulate over
// - Set detector type (0 - no detector, so particles are earth emerging; 1 - spherical detector, particles saved upon entering volume
//    2 - cylindrical detector)
// - Bool to decide if neutrinos and leptons are saved
// - Set threshold energy
//-----------------------------------------------------------------------------
// calling the code is the form ./Simu_elost 1E+20 95.0 1E+4 0 0 4.0 0.92 15
// Command line parameters
// ./Simu_elost <Initial Neutrino energy> <Emergence Angle> <Number of neutrinos> <Cross section model> <Cross section model> 
//              <Energy loss model> <ice depth> <ice density> <particle type>
// Particle types follow pythia particle codes
// +=Particle, -=Anti Particle. 15=tau, 16=nutau, 13=muons, 14=numus, 12=nue,11=e


//updated on 2/18/22

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


using namespace std;


/*
typedef struct {
  int	 nevt; 	      // Event number
  int  ncc; 	      // number of CC interaction suffered ?
  int  nnc;	        // number of NC interaction suffered ?
  int  ndk;	        // number of tau decays
  int  npart;       // Number of particles created (tau or nu_tau only)
  int  trig;        // trigger = 1 if tau finally emerges from Earth
  int  id[40];      // id of produced particle: id=0 if tau neutrino, id=1 if tau
  double theta;     // zenithal angle of initial nu_tau
  double Lmax1;     // Earth thickness crossed
  double Lmax2;     // Earth thickness crossed
  double L0[40];    // Point of interaction of initial neutrino ??
  double Estart;    // Neutrino energy at start point of propagation
  double Eend;      // Particle energy (neutrino or tau) that emerges from Earth
  double E1[40];    // For each particle created: energy at creation
  double E2[40];    // For each particle created: energy of decay (if tau), energy at interaction (if nu)
  double v1[40];    // For each particle created: depth of creation
  double v2[40];    // For each particle created: depth of disappearance
                    // Note: a tau neutrino created in a NC interaction is regarded as a new particle
  double Shheight;  // For a tau emerging, distance traveled before decaying
  double Shlong;
} MYEVT_DEF;
*/
typedef struct {
  stack<int> parent_part_ind;
  stack<int> part_type;       // pythia code - only positive values
  stack<double> part_energy;  // E_min<part_energy<E_init
  stack<double> part_pos;     // position inside earth
  stack<int> anti;           // true if anti particle
  stack<int> generation;      // number of interactions before this particle was produced
  stack<int> NC_num;          // number of NC preceding creation of this particle
  stack<int> CC_num;          // number of CC preceding creation of this particle
  stack<int> GR_num;
  stack<int> dc_num;          // number of decays it took to get to the current particle
  stack<double> x_pos;
  stack<double> y_pos;
  stack<double> z_pos;
  stack<double> traversed_gram;

} particle_info_def;

typedef struct {

  int tau_type[100000][6];       // hold particle types created in tau decay
  double tau_energy[100000][6];  // hold energy of created particle  in tau decay
  int mu_type[100000][6];        // hold particle types created in muon decay
  double mu_energy[100000][6];   // hold energy of created particles in muon decay
    
} reaction_tables_def;

typedef struct {
  string data_dir;            //directory to store the data files
  int anti;                   //+1 is normal, -1 if anti particle
  int starting_type;          // starting type of neutrinos
  bool regen;                 // true or false for lepton regeneration through decays
  bool conversion;            // true or false to include produced particles from W +/- decays
  bool energy_distribution;   // true or false to use energy distribution
  int detector;               // decides detector type, 0=Earth Emerging, 1=Spherical Volume, 2=Cylindrical Volume
  double det_volume;          //Det volume used to calculte radius for spherical model
  bool save_neutrinos;        //Choose is neutrinos are saved
  bool save_charged;          //Choose is charged leptons are saved
  double energy_threshold;    //Set energy threshold
  bool save_events;           //choose to save events
  

}config_init;   // data struct to hold the value read from config file

typedef struct
{
  
  double xi;
  double yi;
  double zi;
  double xf;
  double yf;
  double zf;
  double xe;
  double ye;
  double ze;
  double eex;
  double eey;
  double eez;
  //int input_types;
  //double input_e;
}input_file;//input file parameters setting startpoint, end point, and etrance point to the volume


typedef struct 
{
  double inner_rad;
  double outer_rad;
  double inner_depth;
  double outer_depth;

}det_geom;


//MYEVT_DEF event;
particle_info_def particle_data;
reaction_tables_def reaction_data;
config_init config;
det_geom det;

bool in_volume(double x,double y,double z);


string make_particle_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp);
string make_event_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp);
// initialize reaction from pythia table and convert pythia type to tag code used in code
void initialize_reaction(int tau_type[][6], int mu_type[][6],double tau_ene[][6],double mu_ene[][6]);

// converts pythia tags to tags in this code. converts anti particles to normal matter
int convert_types(int pythia_type);

//load values from config file
void load_config();
void load_geo();
int load_input(input_file *in,int length, string filename);
// -------------------------------------------------
// For lepton energy loss: dE/dX = -alpha + beta(E)*E 
double funcalph(double *x, int *par, int type);

double delta(double X);

// Parameterisations for beta 
double beta9fit(double *x, int *par, int ELOSSmode, int type);	

// Elost by tau and muon dE/dX in GeV/(g/cm^2)
double elost(double E, double dens, int ELOSSmode, int type);

// -------------------------------------------------
// Probability of tau and muon lepton decay
double dPdesdx(double E,int type);	
	
// -------------------------------------------------
// Tau and muon neutrino cross sections: CC,NC,GR
double dsigCC(double E, int CCmode, int type, int AntiNu);	
double dsigNC(double E, int CCmode, int type, int AntiNu);
double dsigGR(double E, int type, int AntiNu);

// -------------------------------------------------
// Local density as a function of zenith angle
double earthdens( double *x, double *par); 
//Local density as a function of x,y coordinates
double get_dens_from_coords(double *coords);

// -------------------------------------------------
// Average density as a function of zenith angle
double mean( double *x, double *par);
double mean_dens_chord(double theta);

//---------------------------------------------------------------
void make_dirs(string dirs);
// Initialize Earth class
// The arguments are water thickness and density. 
// They are initialized to bare rock here but it is re-initialized below.
Earth *terra = new Earth(0.0, 2.6); 

//#############################################################
// Main code
//#############################################################
int main(int argc, char **argv)
{
  double time_start=time(NULL); //set timing variable
  
  
  // call function to load config
  load_config(); 
  make_dirs(config.data_dir);
  load_geo();//might not be necessary now that points are defgined outside of sim
  

  string in_file="in_files/Example_Trajectories (4).csv";
  int input_num = 10003;//will need to change so we can programmatically find the number

  input_file input[input_num];
  
  load_input(input,input_num,in_file);
  
  
  // initializes arrays to hold decay products and populates them from pythia file
  for(int i=0;i<100000;i++){for(int j=0;j<6;j++){reaction_data.tau_type[i][j]=reaction_data.mu_type[i][j]=-1;reaction_data.tau_energy[i][j]=reaction_data.mu_energy[i][j]=0.0;};}  
  initialize_reaction(reaction_data.tau_type,reaction_data.mu_type,reaction_data.tau_energy,reaction_data.mu_energy);
  double initial_energies[100000];
  double neu_energies[100000];
  double lep_energies[100000];
  int gr_counter=0;
  //cout<<config.detector<<" "<<config.det_volume<<endl;
 

  int type_to_save[5]={-1,-1,-1,-1,-1}; //max of 5 types since electrons are forgotten
 
  if(config.save_charged){type_to_save[0]=13;type_to_save[1]=15;} //add muons and taus to pareticle type to savce
  if(config.save_neutrinos){type_to_save[2]=12;type_to_save[3]=14;type_to_save[4]=16;} //add neutrinos to type to save
  
  //for(int i=0;i<5;i++) cout<<type_to_save[i]<<" "; //print which particles will be saved
  //cout<<endl;

  cout << "Lepton Propagation code" << endl;
  // ARGUMENTS:
  // argv[1]:   energy in eV (e.g. 1.e20)
  // argv[2]:   angle (e.g. 91.7) NOTE: This is 180 - exit_angle
  // argv[3]:   (optional) number of events (default:  1e7)
  
  // argv[4]:  (optional) cross-section mode (default: 0 middle, other options: 1 lower, 2 upper.)
  // argv[5]:  (optional) Eloss mode (default: 0 ALLM, other options: 1 ASW)
  
  // argv[6]:  (optional) Water layer thickness in km.   Default value is 0.0 km
  // argv[7]:  (optional) Water layer density in g/cm^3. Default value is 2.6 g/cm3 (rock)
  
  // argv[8]:  (optional)  output file tag (e.g. 0p0km_ice)
  // argv[9]:  (optional)   outputs directory
  // argv[10]: (optional) the directory this program lives in (only needed for cluster runs)
  
  
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
    //tau above and muon below
    (void)strcpy(muondata, "tables/nu_mu_samples");
    (void)strcpy(mfinalccfile, "tables/final_cteq5_cc_nu.data");
    (void)strcpy(mfinalncfile, "tables/final_cteq5_nc_nu.data");
    (void)strcpy(mfinalccbarfile, "tables/final_cteq5_cc_nubar.data");
    (void)strcpy(mfinalncbarfile, "tables/final_cteq5_nc_nubar.data");
  }
  
  //if(argc>10){
  //  (void)strcpy(taudata,     argv[10]);
  //  (void)strcpy(finalccfile, argv[10]);
  //  (void)strcpy(finalncfile, argv[10]);
  //  (void)strcat(taudata, "/Tables/tau_decay_tauola.data");
  //  (void)strcat(finalccfile, "/Tables/final_cteq5_cc_nu.data");
  //  (void)strcat(finalncfile, "/Tables/final_cteq5_nc_nu.data");
  //}
  /*
  cout << "\nTable File Names:" << endl;
  cout << "taudata     " << taudata << endl;
  cout << "tfinalccfile " << tfinalccfile << endl;
  cout << "tfinalncfile " << tfinalncfile << endl << endl;

  cout << "\nTable File Names:" << endl;
  cout << "muondata     " << muondata << endl;
  cout << "finalccfile " << mfinalccfile << endl;
  cout << "finalncfile " << mfinalncfile << endl;
  cout << "finalccbarfile " << mfinalccbarfile << endl;
  cout << "finalncbarfile " << mfinalncbarfile << endl << endl;
  */
  int InitTau = TauData.InitTable(taudata);
  FinalTable *tCCFinalData = new FinalTable;
  FinalTable *tNCFinalData = new FinalTable;
  
  int InitMuon = MuonData.InitTable(muondata);
  FinalTable *mCCFinalData = new FinalTable;
  FinalTable *mNCFinalData = new FinalTable;
  FinalTable *mCCBarFinalData = new FinalTable;
  FinalTable *mNCBarFinalData = new FinalTable;
  
  tCCFinalData->InitTable(tfinalccfile);
  tNCFinalData->InitTable(tfinalncfile);
  mCCFinalData->InitTable(mfinalccfile);
  mNCFinalData->InitTable(mfinalncfile);
  mCCBarFinalData->InitTable(mfinalccbarfile);
  mNCBarFinalData->InitTable(mfinalncbarfile);  

  //bool useEnergyDistribution = false;
  if (atof(argv[1]) == 0) {
    config.energy_distribution = true;
    cout << "Will throw uniformly random x neutrinos energy between log10(E_nu/eV) = 15 and  log21(E_nu/eV)" << endl;
  }

  
  // Initialize Random number generator.
  struct timeval time_struct;
  gettimeofday(&time_struct,NULL);
  srand((time_struct.tv_sec * 1000) + (time_struct.tv_usec / 1000));

  //cout << "Check for random number uniqueness." << endl;
  //for (int ii = 0; ii<3; ii++)
  //{
  //    cout << "Random Test " << ((double) rand() / (double)(RAND_MAX)) << endl;
  //}
  
 

  double angle_time_start=time(NULL);
  double angle=atof(argv[2]);
  if(config.detector==0 && angle<90)
  {
    cout<<"down going events not available for point exit"<<endl;
    return 0;
  }
  
  int produced_muons=0;
  // Cut in energy below which particles are no longer propagated
  double Elim_eV=config.energy_threshold;       // eV
  double Elim=Elim_eV*1.e-9;  // GeV
  
  //-------------------------------------------------
  // Declare several variables
  double rndm;           // random number throughout the code
  double depth=0; //center of detector
  double Lmax;
  double dL=0.;	 	  // Propagation step along chord length. Initialized to zero but it varies for each step.
  double frac=1e-3;	  // Fraction of energy lost by tau in each step. This value stays constant
  double dPdes;         // Probability of tau decay in dL
  double traversed_grammage; // Track the grammage traversed upon entering the Earth.
  double Energy_GeV;	  // Particle energy
  double Bjorken_y;     // Bjorken y value for interactions CC & NC
 

 //-------------------------------------------------//
 // set geomtries for different detectors
 /*
  double det_entrance_d=0;
  double det_exit_d=0;
  double small_r=0;
  if(config.detector==1) //sphere
  {
    depth=2e5;
    small_r=cbrt(3*config.det_volume/4/3.14159);
    small_r=small_r*pow(10,5); //convert km to cm
    det_entrance_d=-small_r;
    det_exit_d=small_r;
    cout<<"here?"<<endl;
  }
  double radius=15*pow(10,5); //ARA
  double height=2.8*pow(10,5); //ARA
  double trans_angle=180/PI*atan(height/2/radius);
  double dh=0;
  double dr=0;
  if(config.detector==2) //cylinder
  {
    depth=2e5;
    if(angle<trans_angle ||angle>(180-trans_angle)) 
    {
      dr=abs(height/2*tan(PI/180*(180-angle)));
      small_r=sqrt(radius*radius+dh*dh);
    }
    if(angle>=trans_angle &&angle<=(180-trans_angle))
    {
      dh=abs(radius*tan(PI/180*(angle-90)));
      small_r=sqrt(height/2*height/2+dr*dr);
    }
    cout<<"or here?"<<endl;
  }
  if(depth<small_r) //just a check for detector limites
  {
    cout<<"ending code, detector would be above surface"<<endl;
    return -1;
  }

  //-------------------------------------------------//
  //set geometry variables for simulation

  //cout<<"small r is "<<small_r<<endl;
  double x_step=cos((angle-90)*PI/180);
  double y_step=sin((angle-90)*PI/180);
  //set start and end points *now capable of looking at down going events too!!
  double start_point[2];
  double end_point[2];
  
  end_point[0]=small_r*x_step;//puts end point on left side
  end_point[1]=R0-depth+small_r*y_step;
  cout<<"depth is "<<depth<<endl;
  //finding start point
  double point_slope=0;
  double y_pos=0;
  double y_neg=0;

  //set for specific cases where tan function causes issues
  if(angle!=180.0 || angle!=0.0) 
  {
    point_slope=-tan((angle-90)*PI/180);
    double det=4*(R0-depth)*(R0-depth)-4*(1+point_slope*point_slope)*((R0-depth)*(R0-depth)-R02*point_slope*point_slope);
    double denom=2*(1+point_slope*point_slope);
    if(det<0)
    {
      cout<<"invalid, breaking code"<<endl;
      return -1;

    }
    y_pos=(2*(R0-depth)+sqrt(det))/denom;
    y_neg=(2*(R0-depth)-sqrt(det))/denom;
    //cout<<"y pos is "<<y_pos<<" , "<<"y neg is "<<y_neg<<endl;
  }
  if(angle>=90.0) start_point[1]=y_neg;
  if(angle<90.0) start_point[1]=y_pos;
 
  
  start_point[0]=sqrt(R02-start_point[1]*start_point[1]);
  if(angle==180.0)
  {
    start_point[1]=-R0;
    start_point[0]=0;
  }
  if(angle==0)
  {
    start_point[1]=R0;
    start_point[0]=0;
  }
  cout<<"start point is "<<start_point[0]<<" , "<<start_point[1]<<endl;
  cout<<"end point is "<<end_point[0]<<" , "<<end_point[1]<<endl;
  
  //sets new max distance for loop... Is set to end right at the beginning of the volume or surface
  */
  
  //cout<<maxL<<" cm"<<endl;
  


  //int brkcnt=0, brkcnt2=0; // These are used to flag particles that have fallen below the energy threshold (Elim)
  //int prop_mode=0;	  // Used in tau propagation
  


  double finalstatecc[2],finalstatenc[2]; // will contain Bjorken (1-y) and y
  
  float dens; // local density along trajectory
  

  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  // Change following parameters depending on application
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  // Number of neutrinos simulated at each zenith angle
  int tot_evt=1e5;
  if(argc>3){
    tot_evt=(int)atof(argv[3]);
  }
  /*
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  cout << "======================================" << endl;
  cout << "Number of neutrinos simulated = " << tot_evt << endl;
  cout << "Type of initial neutrinos = "<<config.starting_type<<endl;
  cout << "Energy = " << atof(argv[1]) << " eV" << endl;
  cout << "Threshold energy = " << Elim*1.e9 << " eV" << endl;
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  */

  //for rounding energies and angle
  string e_temp=to_string(log10((double)atof(argv[1])));
  string es_temp="";
  string ang_temp=to_string(angle);
  string angs_temp="";
  int count=0;
  for(int i =0;i<4;i++)es_temp+=e_temp[i];
  
  if(angle>=100) count=6;
  else count =5;
  for(int i=0;i<count;i++) angs_temp+=ang_temp[i];
    

  

  //cout<<angs_temp<<endl;
  //-------------------------------------------------
  // Output file names using input arguments
  string nameEnergies="";
  string nameEvents="";
  nameEnergies=make_particle_dir(argc,argv,config.data_dir,es_temp,angs_temp); 
  nameEvents=make_event_dir(argc,argv,config.data_dir,es_temp,angs_temp); 
  string nameInFlux="ARA/particles/InFlux_"+es_temp+".dat";
  string nameOutFlux="ARA/particles/OutFlux_"+es_temp+".dat";
  /*
  if(argc>=9){
    nameEnergies+=argv[9];
  }
  */
  //name output file for particles
  
  
  /*
  nameEnergies+=config.data_dir;
  nameEnergies+="/particles";
  nameEnergies+="/particles_";
  nameEnergies+=es_temp;
  nameEnergies+="_";
  nameEnergies+=angs_temp;
  
  if(argc==8){
    nameEnergies+="_";

    nameEnergies+=argv[6];
    nameEnergies+="km_ice_";
    if(atof(argv[4])==0) nameEnergies+="mid";
    if(atof(argv[4])==1) nameEnergies+="low";
    if(atof(argv[4])==2) nameEnergies+="upp";
    nameEnergies+="CS_";
    if(atof(argv[5])==0) nameEnergies+="std";
    if(atof(argv[5])==1) nameEnergies+="low";
    nameEnergies+="EL";
  };
  nameEnergies+=".dat";
  if(argc>8) config.starting_type=atoi(argv[8]);
  */
  ofstream outEnergies(nameEnergies.c_str());
  outEnergies << "10^"<< log10(tot_evt)<<" initial neutrinos of type "<<config.starting_type<<" at energy "<<argv[1]<<". \ntype, anti, NC, CC, GR, DC, Gen, InitNuNum, InitNeutrinoType, OutEnergy, InitEnergy, Part_Pos.\n";
  //cout<<nameEnergies<<endl;
  ofstream outEvents(nameEvents.c_str());
  outEvents<<"vert_x,vert_y,vert_z,tra_x,tra_y,tra_z,E_nu,inel,part_type,i_type,traj_num,p_thrown"<<endl;

  ofstream inFlux(nameInFlux.c_str());
  ofstream outFlux(nameOutFlux.c_str());
  inFlux<<"vert_x,vert_y,vert_z,tra_x,tra_y,tra_z,part_type,anti,energy,part_num,traj_num"<<endl;
  outFlux<<"vert_x,vert_y,vert_z,tra_x,tra_y,tra_z,part_type,anti,energy,part_num,traj_num"<<endl;

  // Get cross-section mode to use
  int CCmode = 0;
  if(argc>4) CCmode = atoi(argv[4]);
  
  // Get cross-section mode to use
  int ELOSSmode = 0;
  if(argc>5) ELOSSmode = atoi(argv[5]);
  
  // Set water layer properties (bare rock is default)
  if(argc>6){
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
  double refTheta=angle;
  //cout << "Theta " << refTheta << " deg" << endl;
  
  //double min_depth = R0-R0*sin(PI*(1.-refTheta/180.));
  //uncomment for icecube
  /*
  if(Depth2>min_depth) 
  {
    cout << "Angle too steep to be measured by IceCube" << endl;
    break;
  }
  */
  
  Lmax=2.*R0*cos(PI*(1.-refTheta/180.));
  
  // Average Earth density
  cout << "Average Earth density " << mean_dens_chord(refTheta) << " g/cm^3" << endl;
  cout << endl;
  
  cout << "======================================" << endl;
  cout << "Emerging Lepton Energy file: " << endl;
  cout << nameEnergies << endl;
  cout << endl;
  
  cout << "======================================" << endl;
  

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //                                    Initiate run
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //double maxL=sqrt((start_point[0]-end_point[0])*(start_point[0]-end_point[0])+(start_point[1]-end_point[1])*(start_point[1]-end_point[1]));
  //if(Lmax2<=0.) Lmax2=0.;                 // set to zero  if the value is negative.
  if(Lmax<=0.) Lmax=0.;  
  //printf ("Lmax %1.3e cm, %1.2f deg \n", Lmax2, refTheta);//change to lmax2 for icecube
  
  // create look-up datable of distance as a function of grammage
  
  // start with one particle then it will loop over stack and iterate to the next original particle
  int new_muons=0; //count of pair of particles made. +1 for muon +1 for neutirno
  int total_CC=0;
  int total_NC=0;
  int total_GR=0;
  double num_tau_decays=0;
  double num_muon_decays=0;
  double distance_loop_num=0;
  int tau_neutrinos_below=0;
  int muon_neutrinos_below=0;
  int muons_below=0;
  int taus_below=0;
  bool split_type=false;
  int part_count=0;
  if(config.starting_type==0) split_type=true;
  //cout<<"help"<<endl;
  //cout<<small_r<<", "<<Lmax<<endl;

  //start loop over input trajectories
  int event_count=0;
  bool run_multi=true;
  int sim_count=0;
  
  while(sim_count<10 && run_multi)
  {
    sim_count++;
    run_multi=false;
  for( int i=0;i<input_num;i++)//loop over trajectories
  { 

    int num_count=0;
    //---------------------------------
    double sum_grammage = 0.;
    
    double check_grammage=0;
    double maxL=sqrt((input[i].eex-input[i].xi)*(input[i].eex-input[i].xi)
                    +(input[i].eey-input[i].yi)*(input[i].eey-input[i].yi)
                    +(input[i].eez-input[i].zi)*(input[i].eez-input[i].zi));
    if(maxL==0)
    {
      cout<<"didnt load correctly";
    }
    double x_step=(input[i].eex-input[i].xi)/maxL;
    double y_step=(input[i].eey-input[i].yi)/maxL;
    double z_step=(input[i].eez-input[i].zi)/maxL;
    
    double temp_pos[3]={input[i].xi,input[i].yi,input[i].zi};
    
    for (int ii=1; ii<=1000000; ii++) //reduce by 10
    {
    double dl = maxL/1000000; //change to Lmax2 for icecube
    double x_val = maxL * double(ii) / 1000000.;//change to lamx2 for icecube
    
    temp_pos[0]=temp_pos[0]+dl*x_step;
    temp_pos[1]=temp_pos[1]+dl*y_step;
    temp_pos[2]=temp_pos[2]+dl*z_step;
    sum_grammage +=dl*get_dens_from_coords(temp_pos);
    //sum_grammage +=  dx*earthdens(&x_val, &Lmax);//change to lamx2 for icecube
    
    }
    /*
    for (int ii=1; ii<=1000000; ii++)
    {
    double dx = maxL/1000000;
    double x_val = maxL * double(ii) / 1000000.;
    check_grammage +=  dx*earthdens(&x_val, &maxL);
    }
    cout<<"new sum grammage "<<sum_grammage<<". old sum grammage" <<check_grammage<<endl;
    return 0;
    //printf("sum_grammage %1.5e g/cm^2\n",sum_grammage);
    //printf("mean_dens %1.2e\n\n",sum_grammage/Lmax2);//chnage to lmax2 for icecube
    */
    double* cumulative_grammage = new double[1000000];
    double* grammage_distance = new double[1000000]; // in cm
    
    double d_grammage = sum_grammage/1000000.; // g/cm^2
    
    cumulative_grammage[0] = 0.;
    grammage_distance[0]   = 0.;
    temp_pos[0]=input[i].xi;
    temp_pos[1]=input[i].yi;
    temp_pos[2]=input[i].zi;


    for (int ii=1; ii<=1000000; ii++)
    {
    
    double dl = d_grammage/get_dens_from_coords(temp_pos);//chnage to lmax2 for icecube
    double l_val = grammage_distance[ii-1];
    cumulative_grammage[ii] = cumulative_grammage[ii-1] + dl*get_dens_from_coords(temp_pos);//change to lamx2 for icecube
    grammage_distance[ii] = l_val + dl;  
    temp_pos[0]=temp_pos[0]-dl*x_step;
    temp_pos[1]=temp_pos[1]+dl*y_step;
    temp_pos[2]=temp_pos[2]+dl*z_step;
    //printf("*** ii %d %1.5f %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
    //if(ii%100000 ==0) printf("ii %d %1.2e %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
    }


    //save charged leptons entering the volume!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    //----------------------------------
    cout<<"primary trajectory: "<<i<<endl;
    bool got_event=false;//!got_event && num_count<1000 line below
    while(num_count<1000)//!got_event for running until getting 1 event or num_count<N for throwing N particles
    {
    num_count++;
    
    double pos[3]={input[i].xi,input[i].yi,input[i].zi};
    double en_pos[3]={input[i].xe,input[i].ye,input[i].ze};
    double len_ent=sqrt((input[i].xe-input[i].xi)*(input[i].xe-input[i].xi)
                    +(input[i].ye-input[i].yi)*(input[i].ye-input[i].yi)
                    +(input[i].ze-input[i].zi)*(input[i].ze-input[i].zi));
    bool has_been_saved=false;
    bool save_cond=false;
    bool exit_cond=false;
    //cout<<"primary particle "<<num_count<<endl;
    //cout<< len_ent<<","<<maxL<<endl;
    //if(len_ent>maxL)cout<<endl<<endl<<endl;
    //cout<<x_step<<","<<y_step<<","<<z_step<<endl;
    //cout<<input[i].xi<<","<<input[i].yi<<","<<input[i].zi<<endl;
    //cout<<input[i].xf<<","<<input[i].yf<<","<<input[i].zf<<endl;
    Energy_GeV = atof(argv[1])*pow(10,-9); // Get nu_tau energy from input argument
    if (config.energy_distribution) Energy_GeV =  pow(10,6 + (6 * (double)rand()/RAND_MAX));
    //cout<<"initial energy GeV is "<<Energy_GeV<<endl;
    //cout<<"threshold energy in GeV is "<<Elim<<endl;
    int part_type=0;
    
    if(config.starting_type==16) part_type=16;
    if(config.starting_type==14) part_type=14;
    if(config.starting_type==12) part_type=12;

    if(split_type==true && part_count<tot_evt/3) part_type=12;
    if(split_type==true && part_count>=tot_evt/3 && part_count<=2*tot_evt/3 ) part_type=14;
    if(split_type==true && part_count>2*tot_evt/3) part_type=16;
 
    if(part_type==0)
    {
      cout<<"Particle type failed to intialize...breaking code\n";
      return 0;
    }

    part_count++;
    double part_energy=Energy_GeV;
    double part_pos=0;
    int anti=1;
    int generation=0;
    int NC_num=0;
    int CC_num=0;
    int GR_num=0;
    int dc_num=0;
    anti=config.anti;

    //stuff for w decays

    const double e_rest=(mW*mW+me2)/(2*mW);
    const double m_rest=(mW*mW+mmuon2)/(2*mW);
    const double t_rest=(mW*mW+mtau2)/(2*mW);
    
  
    traversed_grammage = 0.; //initialize traversed grammage for each event

    
    int loop_num=0;
    //======================= Loop over stacks until empty
    do 
    { 
      
      loop_num++;
      
      has_been_saved=false;
      bool save_cond=false;
      bool exit_cond=false;

      //pop all particle data from the stacks to working variables
      //======================================================
      if(!particle_data.part_type.empty())
      {
        part_type=particle_data.part_type.top();
        part_energy=particle_data.part_energy.top();
        part_pos=particle_data.part_pos.top();
        anti=particle_data.anti.top();
        generation=particle_data.generation.top();
        NC_num=particle_data.NC_num.top();
        CC_num=particle_data.CC_num.top();
        dc_num=particle_data.dc_num.top();
        GR_num=particle_data.GR_num.top();
        pos[0]=particle_data.x_pos.top();
        pos[1]=particle_data.y_pos.top();
        pos[2]=particle_data.z_pos.top();
        traversed_grammage=particle_data.traversed_gram.top();
        //cout<<part_type<<endl;
        //======================================================
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
      }
      double tau_path=0;
      //if(loop_num!=0)cout<<"stack loop of type "<<part_type<<" and energy "<<part_energy<<"at x,y "<<pos[0]<<","<<pos[1]<<endl;
      if(part_type==11) {continue;} //ignores particles of electron flavor
      if(part_energy<Elim) continue; //ignore particles below threshold in case they make it through
      
      //cout<<part_type<<","<<anti<<endl;

      // Flag to see how fast code is running
      //if(!((float)i/100000-(int)(i/100000))) cout<< i << endl;
      

      //brkcnt=0;
      //prop_mode =0;
      
      //======================= Start propagation along chord of length Lmax
      //printf("Ldist, Lmax %1.2e %1.2e\n", Ldist, Lmax);
      //Ldist=0.;
      
      bool broken=false;
      bool left_volume=false;
      bool entered_volume=false;
      while(part_pos<maxL &&  !left_volume) //EDIT HERE FOR END CONDITION
      {
        
        bool in_vol=in_volume(pos[0],pos[1],pos[2]);
        //create holding arrays for reactions
        int reaction_types[6]={0,0,0,0,0,0};
        double reaction_energies[6]={0,0,0,0,0,0};
        int anti_type[6]={1,1,1,1,1,1};
        //cout <<"distance is "<< part_pos << "and energy is "<<part_energy<<endl;
        // Get the local density for this pa rt of the chord.
        dens = get_dens_from_coords(pos);
        if(dens==0){cout<<"local dens is zero"<<endl;}

        //cout<<"In dist loop - type: "<<part_type<<". energy: "<<part_energy<<". pos: "<<part_pos<<". gen: "<< generation<<endl;
        bool change=false; //is needed so in case a NuTau becomes a Tau it won't be seen by the "is tau" section ie skipping a step

        //===========================
        // Particle is a neutrino
        //===========================
        if(part_type==12||part_type==14||part_type==16) //CHANGE PARTICLE TYPE TO PYTHIA
          {

          // Number of interaction lengths propagated in this step is given by an exponentially distributed random number.
          double num_int_lens=-log((double) rand() / (double)(RAND_MAX)); // Randomly sampled number of interaction lengths.
          
          // Convert the number of interaction lengths to grammage.
          double X_int;
          X_int = num_int_lens /(Navo*(dsigGR(part_energy,part_type,anti)+dsigCC(part_energy, CCmode,part_type,anti)+dsigNC(part_energy, CCmode,part_type, anti)));
          
          // The following lines use the grammage_distance lookup table to estimate what position along the trajectory this Xint corresponds to.
          
          // Add X_int to the total grammage traversed by the particle
          traversed_grammage += X_int;
          
          // Initialize the interaction length distance for this step to zero.
          //double Lint = 0.;
        
          // If too large, make sure it exits the volume.
          // NOTE: use floats for this condition. Using ints is bad if float > 2^32, then you get negative int.
          if( traversed_grammage/d_grammage + 1. > 1000000.){
            part_pos = maxL; // chnage to lamx2 for icecube NOTE: 1000000. is the size of the look-up table.
            traversed_grammage = sum_grammage;
            pos[0]=input[i].xf;
            pos[1]=input[i].yf;
            pos[2]=input[i].zf;
            
          }
          // If contained within the trajectory, linearly interpolate its interaction distance.
      
          if ( floor(traversed_grammage/d_grammage) + 1. < 1000000.) // NOTE: 1000000. is the size of the look-up table.
          {
            double before_step=part_pos;
            // Get the entry in the look-up table corresponding to the traversed grammage
            int ii_grammage = int(traversed_grammage/d_grammage) + 1;
      
            // Linearly interpolate to estimate the distance propagated
            double slope = (grammage_distance[ii_grammage] - grammage_distance[ii_grammage-1])/d_grammage;
          
            double intercept = grammage_distance[ii_grammage] - slope*cumulative_grammage[ii_grammage];
            //Lint = slope*traversed_grammage + intercept - part_pos; // keep track of this step's interaction length.
            part_pos = slope*traversed_grammage + intercept ;
            double after_step=part_pos;
            pos[0]=pos[0]-(after_step-before_step)*x_step;
            pos[1]=pos[1]+(after_step-before_step)*y_step;
            pos[2]=pos[2]+(after_step-before_step)*z_step;
          }
          if(in_volume(pos[0],pos[1],pos[2])&& !entered_volume&&config.save_neutrinos)
          {
            entered_volume=true;
            inFlux<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<x_step<<" "<<y_step<<" "<<z_step<<" "<<part_type<<" "<<anti<<" "<<part_energy<<" "<<num_count<<" "<<i<<endl;
          }
        
          //cout<<889<<",";
          // Save the interaction distance propagated in event structure
          //event.L0[npart]=Lint;
          
          // if the neutrino interaction is still inside Earth, simulate an NC or CC interaction and check that particle is still above the tracking energy threshold (Elim)
          if(part_pos<maxL)// for icecube  part_pos<Lmax //!((Lmax1<=part_pos) && (part_pos<=Lmax2))&&part_pos<=Lmax2
          {
            // Check if it is CC or NC interaction
            double total_cross=dsigNC(part_energy, CCmode,part_type,anti)+dsigCC(part_energy, CCmode,part_type,anti)+dsigGR(part_energy,part_type,anti);
            double rand_val=((double) rand() / (double)(RAND_MAX));
            if(total_cross==0){cout<<"total_cross is zero"<<endl;}

            double CCratio=dsigCC(part_energy, CCmode,part_type,anti)/total_cross;
            double NCratio=dsigNC(part_energy, CCmode,part_type,anti)/total_cross;
            double GRratio=dsigGR(part_energy,part_type,anti)/total_cross;
            //if(part_type==12&&anti==-1)cout<<GRratio<<endl;
            bool CChappens=rand_val<CCratio;
            bool NChappens=(CCratio<=rand_val && rand_val<CCratio+NCratio);
            bool GRhappens=(CCratio+NCratio<=rand_val && rand_val<CCratio+NCratio+GRratio);
            

            if(CChappens)
            {
              //=======================
              // CC interaction occurs (the tracked particle changes from tau neutrino to tau lepton.)
              //=======================
              //total_CC++;

              CC_num++;
              //cout<<"charged current of neutrinos happened \n";
              // Obtain Bjorken y
              if(part_type==16)tCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
    
              if(part_type==14)
              {
                if(anti==1)mCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                if(anti==-1)mCCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
              }
              if(part_type==12)
              {
                broken=true;
                part_type=12;
                //particle becomes e so save in event and break
              }

              Bjorken_y=finalstatecc[1];//INEASTIVITY?
            
              // Set the tau lepton energy from the sampled Bjorken y.
              double initial_energy=part_energy;
              part_energy=(1.-Bjorken_y)*part_energy;
              double shower_energy=initial_energy-part_energy;
              // Increment the cc interaction counter in the event structure.
              //event.ncc++;
              if(in_volume(pos[0],pos[1],pos[2]))
              {
              got_event=true;
              event_count++;
              outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_energy<<","<<Bjorken_y<<","<<part_type*anti<<","<<0<<","<<i<<","<<num_count<<endl;
              }
              if(broken==true)break;
              // Increment the particle counter in the event structure
              //event.npart++;
              //npart++;
              generation++;
              // Save the new particle energy and distance in the event structure/
              //event.E1[npart]=Energy_GeV;
              //event.v1[npart]=Ldist;
              
              // Change the particle tag from tau neutrino to lepton.
              //event.id[npart]=1;
              //tag=1;
              
             
              if(part_type==14) part_type=13;//nu_mu -> mu
              if(part_type==16) part_type=15;//nu_tau -> tau
            

              // get the density at the current location before jumping to the tau lepton part of the loop
              //dens = get_dens_from_coords(pos);//change to Lmax2 dfor icecube
              change=true;//signify that particle type changed so simulation of decay delayed until next distance loop iteration
              //ADD EVENT SAVING!!
            }
            else if (NChappens)
            {
              //=======================
              // NC interaction occurs (the tracked particle remains a tau neutrino with reduced energy.)
              //=======================
              //total_NC++;
              //cout<<"neutral current of neutrinos happened \n";
              // Obtain Bjorken y

              if(part_type==16)tNCFinalData->ThrowFinal(log10(part_energy),finalstatenc);  
      
              if(part_type==14 || part_type==12)
              {
                if(anti==1) mNCFinalData->ThrowFinal(log10(part_energy),finalstatenc);
                if(anti==-1) mNCBarFinalData->ThrowFinal(log10(part_energy),finalstatenc);
              }
              
              Bjorken_y=finalstatenc[1]; //INELASTICITY
              
              // Set the neutrino energy from the sampled Bjorken y.
              double initial_energy=part_energy;
              part_energy=(1.-Bjorken_y)*part_energy;
              double shower_energy=initial_energy-part_energy;
              // Increment the cc interaction counter in the event structure.
              //event.nnc++;           // count NC interaction
              NC_num++;
                                    // Increment the particle counter in the event structure
              //event.npart++;
              //npart++;
              generation++;
              if(in_volume(pos[0],pos[1],pos[2]))
              {
              outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_energy<<","<<Bjorken_y<<","<<part_type*anti<<","<<1<<","<<i<<","<<num_count<<endl;
              event_count++;
              got_event=true;
              }
              // Save the new particle energy and distance in the event structure/
              //event.E1[npart]=Energy_GeV;
              //event.v1[npart]=Ldist;
              
              // Keep the particle tag as neutrino. so no change in part_tupe should happen
              //if(tag==0){event.id[npart]=0;}
              //if(tag==2){event.id[npart]=2;}
              change=true;
              //ADD SAVE EVENT!!
            }
            else if (GRhappens)
            {
              //total_GR++;
              GR_num++;
              double initial_energy=part_energy;
              
              if(part_type!=12 || anti!=-1) 
              {
                cout<<CCratio<<" "<<NCratio<<" "<<GRratio<<endl;
                cout<<"CC happens = "<<CChappens<<" NChappens = "<<NChappens<<" GR happens = "<<GRhappens<<endl;
                cout<<"GR without a valid particle"<<endl;
                continue;
              }
              
              double react_rand=((double) rand() / (double)(RAND_MAX));
              //cout<<react_rand<<endl;
              if(react_rand<0.676)
              { 
                
                //cout<<"W+ decayed to quarks"<<endl;
                if(in_volume(pos[0],pos[1],pos[2]))
                {
                
                outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_energy<<","<<1<<","<<part_type*anti<<","<<2<<","<<i<<","<<num_count<<endl;
                event_count++;
                got_event=true;
                }
                broken=true;
                break;
              }
              else
              {
                
                //cout<<"W+ decayed to laptons"<<endl;
                double lep_rand=((double) rand() / (double)(RAND_MAX));
                double initial_E=part_energy;
                double lepton_mass=0;
                double lepton_type=0;
                int lepton_anti=1;
                
                if(lep_rand<(1./3.))    //change particle types
                {//e made
                  part_type=12;
                  lepton_mass=me;
                  lepton_type=11;
                  
                }
                else if(lep_rand<(2./3.) &&lep_rand>=(1./3.))
                {//mu made
                  part_type=14;
                  lepton_mass=mmuon;
                  lepton_type=13;
                  
                }
                else if(lep_rand>=(2./3.))
                {//tau made
                  part_type=16;
                  lepton_mass=mtau;
                  lepton_type=15;
                  
                }
                //cout<<"GR happened at "<<pos[0]<<" "<<pos[1]<<" and decayed to "<<part_type<<" and "<<lepton_type<<endl;
                //cout<<part_type<<","<<lepton_type<<endl;
                part_energy=((double) rand() / (double)(RAND_MAX))*part_energy*(mW*mW-lepton_mass*lepton_mass)/(mW*mW);
                //gamma=E_nu/mW, E_e_rest=mW^2+ml^2)/(2mW), 
                double angle_rand=(double)rand()/(double)RAND_MAX;
                double lepton_energy=initial_E/mW*(mW*mW+lepton_mass*lepton_mass)/(2*mW)*cbrt(8-8*angle_rand);

                part_energy=initial_E-lepton_energy;
                double gr_inel=0;
                if(part_type*anti==12) gr_inel=(initial_E-part_energy)/initial_E;
                if(in_volume(pos[0],pos[1],pos[2]))
                {
                outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_E<<","<<gr_inel<<","<<part_type*anti<<","<<3<<","<<i<<","<<num_count<<endl;
                event_count++;
                got_event=true;
                }
                //if(gr_counter<100000 && lepton_type==15){
                //  initial_energies[gr_counter]=initial_E;
                //  neu_energies[gr_counter]=part_energy;
                //  lep_energies[gr_counter]=lepton_energy;
                //  gr_counter++;
                //}
                

                if(lepton_type!=11) //throw extra lepton into stack and continue on with the produced nu
                {
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
                }
                
              }
              change=true;//signify that particle type changed so simulation of decay delayed until next distance loop iteration
              
            }
          }// end of 'if(Ldist<Lmax)'  (particle still inside Earth)

        }// end of 'if (tag == 0,1,2)' (particle is a nu_e,nu_m,nu_t)
        
        //=========================
        // Particle is a tau lepton or muon lepton
        //=========================
        if((part_type==11||part_type==13||part_type==15)&&change==false&&broken==false) //change particle types
        {
          // Estimate step length based on Energy, dE/dx, local density, and fraction.
          dL=(part_energy/(dens*elost(part_energy, dens, ELOSSmode,part_type)))*frac;
          //cout << "Ldist " << 1.e-5*Ldist << "  R " << 1.e-5*sqrt(R02 - (Ldist*Lmax)+Ldist*Ldist) << " dens " << dens << endl;
          //if(tau_path==0)cout<<"tau energy start "<<part_energy<<endl;
          tau_path+=dL;
          //cout<<dL<<endl;
          // Check if tau leaves the Earth after dL. If it does then adjust last step
          if(part_pos+dL > maxL) 
          {
            dL=maxL-part_pos;//change tolmax1 for icecube
            //cout<<"tau left before decaying "<<endl;
          }
          // Calculate the traversed grammage
          traversed_grammage+=dL*dens;
          
          // Calculate the probability that the tau lepton will decay along the step dL.
          dPdes=1.-exp(-dL*dPdesdx(part_energy,part_type));
          
          rndm=((double) rand() / (double)(RAND_MAX));
          
          //update particle position
          //if(rndm<dPdes)cout<<"tau to decay here "<<part_pos<<endl;
          part_pos=part_pos+dL;
          pos[0]=pos[0]+dL*x_step;
          pos[1]=pos[1]+dL*y_step;
          pos[2]=pos[2]+dL*z_step;
          //if(rndm<dPdes)cout<<"tau decays here "<<part_pos<<endl;
          if(in_volume(pos[0],pos[1],pos[2])&& !entered_volume&&config.save_charged)
          {
            entered_volume=true;
            inFlux<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<x_step<<" "<<y_step<<" "<<z_step<<" "<<part_type<<" "<<anti<<" "<<part_energy<<" "<<num_count<<" "<<i<<endl;
          }
          // Calculate a random number used for determining whether the tau decays or not in this step.
          
         

          part_energy = part_energy-dL*dens*elost(part_energy, dens, ELOSSmode,part_type);
          // Determine whether the lepton decays or not.
          if(rndm > dPdes)
          {
            //==============================
            // The tau lepton does NOT decay
            //=============================
            //cout<<"no decay"<<endl;
            // Advance the propagation distance by the step dL
            
            // Account for the tau lepton energy lost in the step dL
            
          }
          else
          {
            //=======================
            // The lepton decays
            //=======================
            dc_num++;
            //cout<<"tau path length before decaying is "<<tau_path<<endl;
            //cout<<"tau energy at decay "<<part_energy<<endl;
            //cout<<"dL is "<<dL<<endl;
            //cout<<"dpdes is "<<dPdes<<endl;
            //cout<<"rndm is "<<rndm<<endl;
            
            //if(part_type==6) num_tau_decays++;
            //if(part_type==5) num_muon_decays++;

            // Account for the tau lepton energy lost in the step dL
            
            
            // Save the updated particle energy and location of interaction
            //event.E2[npart]=Energy_GeV;
            //event.v2[npart]=Ldist;
          
            // Get the energy of the neutrino produced in the decay
            generation++;
            //if(tag==1) TauData.ThrowFinal(finalstate);
            //if(tag==3) MuonData.ThrowFinal(finalstate);
            //Energy_GeV=finalstate[0]*Energy_GeV;
            double initial_energy=part_energy;
            int reaction_index=(double)rand()/(double)RAND_MAX*100000;
            //cout<<"reaction index is:"<<reaction_index<<"\n\n\n";
            //if tau
            int initial_particle_type=part_type;
            double dec_inel=0;
            int initial_particle=part_type;
            //cout<<"tau decay"<<endl;
            if(part_type==15)
            { 
              //cout<<"tau decay"<<endl;

              if(config.conversion)
              {
                for(int j=1;j<6;j++)
                {
                  if(reaction_data.tau_type[reaction_index][j]!=0) 
                  { 

                    int anti_ness=1;
                    if(reaction_data.tau_type[reaction_index][j]<0) anti_ness=-1;
                    anti_type[j]=anti*anti_ness;
                    reaction_types[j]=abs(reaction_data.tau_type[reaction_index][j]);
                    reaction_energies[j]=part_energy*reaction_data.tau_energy[reaction_index][j];
                    //if(reaction_data.tau_type[reaction_index][j]==-12)cout<<"ve+ made"<<endl;
                    //cout<<reaction_index<<endl;
                    //cout<<"tau decay with type  "<<reaction_types[j]<< "has anti ness of "<<anti_type[j]<<" and reaction energy is "<<reaction_energies[j]<<endl;
                  }
                }
              
              }
              part_energy*=reaction_data.tau_energy[reaction_index][0];
              //cout<<part_energy<<endl;
              //event.id[npart]=0;
              part_type=16;
            }
            if(part_type==13)
            {
              //cout<<"muon decay"<<endl;
              if(config.conversion)
              {
              for(int j=1;j<6;j++)
                {
                  if(reaction_data.mu_type[reaction_index][j]!=0) 
                  { 

                    int anti_ness=1;
                    if(reaction_data.mu_type[reaction_index][j]<0) anti_ness=-1;
                    anti_type[j]=anti*anti_ness;
                    reaction_types[j]=abs(reaction_data.mu_type[reaction_index][j]);
                    reaction_energies[j]=part_energy*reaction_data.mu_energy[reaction_index][j];

                    //cout<<"muon decay with type is "<<reaction_types[j]<<" and reaction energy is "<<reaction_energies[j]<<endl;
                  }
                }
              
              }
              part_energy*=reaction_data.mu_energy[reaction_index][0];
              //event.id[npart]=2;
              part_type=14; 
            }

            
            double frac_energy_dumped=0;
            for(int frac=0;frac<6;frac++)
            {
              if(reaction_data.tau_type[reaction_index][frac]==0 || reaction_data.tau_type[reaction_index][frac]==11)
              frac_energy_dumped+=reaction_data.tau_energy[reaction_index][frac];
            }
            //double lost_energy=0;
            //for(int i=1;i<6;i++)
            // {
            //  if(reaction_types[i]!=4) continue;
            //  lost_energy+=reaction_energies[i];
            //}
            //double shower_energy=initial_energy-part_energy-lost_energy;
            if(in_volume(pos[0],pos[1],pos[2]))
            {
            outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_energy<<","<<frac_energy_dumped<<","<<initial_particle*anti<<","<<4<<","<<i<<","<<num_count<<endl;
            event_count++;
            got_event=true;
            }
            if(!config.regen) { broken=true; break;}
          }
        }// end of 'if (tag == 1,2,3)' (particle is lepton)
          


    

        //time to pop all the extra particles back in the stack to be looped over after
        if(config.conversion)
        {
          for( int j=1;j<6;j++)
          {
              
            if((reaction_energies[j]>Elim)&&(part_pos<Lmax)&&(reaction_types[j]!=0)) //particles inside earth and above threshold are stacked, change to lam2 for icecube
            {
              //cout<<"pushing new particle to stack \n";
              produced_muons+=1;
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
              }
          }
        }
        if(part_energy<Elim){
          //if(part_type==3) tau_neutrinos_below++;
          //else if (part_type==2) muon_neutrinos_below++;
          //else if (part_type==6) taus_below++;
          //else if (part_type==5) muons_below++;
          broken=true;
          break;
          }
        if(in_vol&&!in_volume(pos[0],pos[1],pos[2]))
          {
            //cout<<"left volume"<<endl;
            left_volume=true;
          }
      
      } // ends loop 'for(Ldist=0.;Ldist<Lmax;)'
      
      //=================================================
      // Write energy of emerging tau to output text file
      //=================================================
      //cout<<"save- type: "<<part_type<<" .energy: "<<part_energy<<".pos: "<<part_pos<<". gen: "<< generation<<endl;
      if(left_volume&&config.save_neutrinos)
      {
        
        outFlux<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<x_step<<" "<<y_step<<" "<<z_step<<" "<<part_type<<" "<<anti<<" "<<part_energy<<" "<<num_count<<" "<<i<<endl;
      }
      else if(left_volume&&config.save_charged)
      {
        
        outFlux<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<x_step<<" "<<y_step<<" "<<z_step<<" "<<part_type<<" "<<anti<<" "<<part_energy<<" "<<num_count<<" "<<i<<endl;
      }

      for(int j=0; j<5;j++)
      {
        //cout<<type_to_save[i]<<","<<part_type<<endl;
        //cout<<"part energy is "<<part_energy<<"with threshold "<<Elim<<endl;
        if((part_type==type_to_save[j])&&(part_energy>Elim)&&broken==false&&has_been_saved==false)//changge to just poarticle tyoe if saving all particls even below thrwhpold
          {//add &&false after testing influx and outflux files
            
            outEnergies <<part_type<<" "<<anti<<" "<<NC_num << " " << CC_num << " " << GR_num<<" "<<
            dc_num << " " <<generation <<" "<<i<< " " <<
            setprecision(7)  << log10(part_energy)+9 << " " << log10(Energy_GeV)+9<<" "<<part_pos<<"\n";
            has_been_saved=true;
          }
        
      }//if the main particle isn't saved then it is forgetten falling below thrshold
      
    } while(!particle_data.part_type.empty()); // end of loop of stack
    }
  delete cumulative_grammage;
  delete grammage_distance;
  } // end of loop for initial particles
  

  }
  outEnergies << "END" << endl; // write END in the last line of the text output file. 
  outEnergies.flush();
  outEvents <<"END"<<endl;
  outEvents.flush();
  inFlux<<"END"<<endl;
  outFlux<<"END"<<endl;
  outFlux.flush();
  inFlux.flush();
  double angle_time_end=time(NULL);
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
 
  double time_elapsed=time(NULL)-time_start;
  cout<<"Time for this angle and energy is "<<time_elapsed<<" s"<<endl;
  cout<<event_count<<" events"<<endl;
  cout<<sim_count<<" loops"<<endl;
  return 0;
  
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// ===================================================
// Several functions
// ===================================================
bool in_volume(double x,double y,double z)
{
  double rad=sqrt(x*x+y*y);
  if(rad<=15*pow(10,5)&&(z<R0)&&(z>R0-2.8*pow(10,5)))
  {
    return true;
  }
  else return false;
}
int load_input(input_file *in,int file_length,string filename)
{
  fstream test_csv(filename.c_str());
  string row = "";
  getline(test_csv,row);
  
  double holding[12];
  int row_ind=0;
  double conv=pow(10,5);//km to cm
  while(row.size()!=0)
  {
    
    int i=0;
    int place=0;
    int length=0;
    int index=0;
    while(row.find(",",place)!=-1)
        {
          
          index=row.find(",",place);
          length=index-place;
          holding[i]=conv*atof(row.substr(place,length).c_str());
          place=index+1;
          i++;
          
        }
    holding[11]=conv*atof(row.substr(place,row.size()-place).c_str());
    
    in[row_ind].xi=holding[0];
    in[row_ind].yi=holding[1];
    in[row_ind].zi=holding[2];
    in[row_ind].eex=holding[3];
    in[row_ind].eey=holding[4];
    in[row_ind].eez=holding[5];
    in[row_ind].xf=holding[9];
    in[row_ind].yf=holding[10];
    in[row_ind].zf=holding[11];
    in[row_ind].xe=holding[6];
    in[row_ind].ye=holding[7];
    in[row_ind].ze=holding[8];
    //for(int abs=0;abs<12;abs++)
    //{
    //  cout<<holding[abs]<<endl;
    //}
    //in[row_ind].input_types=(int)holding[9];
    //in[row_ind].input_e=holding[10];
    //from here assign the holding[] to the correct parameter in the sim
    //with a row index ie x[row_ind]=holding[0]
    getline(test_csv,row);
    //cout<<"\n\n";
    row_ind++;
    if(row_ind>file_length) 
    {
      cout<<"reached end of file"<<endl;
      cout<<row_ind<<">"<<file_length<<endl;
      return 0;
    }
  }
  return 0;
}

string make_particle_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp)
{
  string nameEnergies="";
  nameEnergies+=config.data_dir;
  nameEnergies+="/particles";
  nameEnergies+="/particles_";
  nameEnergies+=es_temp;
  nameEnergies+="_";
  nameEnergies+=angs_temp;
  
  if(argc==8){
    nameEnergies+="_";

    nameEnergies+=argv[6];
    nameEnergies+="km_ice_";
    if(atof(argv[4])==0) nameEnergies+="mid";
    if(atof(argv[4])==1) nameEnergies+="low";
    if(atof(argv[4])==2) nameEnergies+="upp";
    nameEnergies+="CS_";
    if(atof(argv[5])==0) nameEnergies+="std";
    if(atof(argv[5])==1) nameEnergies+="low";
    nameEnergies+="EL";
  };
  nameEnergies+=".dat";
  if(argc>8) config.starting_type=atoi(argv[8]);
  return nameEnergies;
}
string make_event_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp)
{
  string nameEvents="";
  nameEvents+=config.data_dir;
  nameEvents+="/events";
  nameEvents+="/events_";
  nameEvents+=es_temp;
  //nameEvents+="_";
  //nameEvents+=angs_temp;
  nameEvents+=".dat";

  return nameEvents;

}


double get_dens_from_coords(double *coords)
{
  double f;
  double radius=sqrt(coords[0]*coords[0]+coords[1]*coords[1]+coords[2]*coords[2]);

  f = terra->GetDensity(radius);

  return f;
}

//loads config.txt file
void load_config()
{
  ifstream fin("config.txt");
  string line;
  while (getline(fin,line))
  {
     istringstream sin(line.substr(line.find("=")+1));
     if((int)line.find("data_dir")!=-1) sin>>config.data_dir;
     else if ((int)line.find("starting_type")!=-1) sin>>config.starting_type;
     else if ((int)line.find("anti")!=-1) sin >>config.anti;
     else if ((int)line.find("regen")!=-1) sin>>config.regen;
     else if ((int)line.find("conversion")!=-1) sin>>config.conversion;
     else if ((int)line.find("energy_distribution")!=-1) sin>>config.energy_distribution;
     else if ((int)line.find("detector")!=-1) sin>>config.detector;
     else if ((int)line.find("det_volume")!=-1) sin >>config.det_volume;
     else if ((int)line.find("save_neutrinos")!=-1) sin >>config.save_neutrinos;
     else if ((int)line.find("save_charged")!=-1) sin >>config.save_charged;
     else if ((int)line.find("energy_threshold")!=-1) sin >>config.energy_threshold;
     else if ((int)line.find("save_events")!=-1) sin >>config.save_events;

     

  }
}

void load_geo()
{
  ifstream fin("det_geom.txt");
  string line;
  while (getline(fin,line))
  {
     istringstream sin(line.substr(line.find("=")+1));
     if ((int)line.find("inner_rad")!=-1) sin>>det.inner_rad;
     else if ((int)line.find("inner_depth")!=-1) sin >>det.inner_depth;
     else if ((int)line.find("outer_rad")!=-1) sin>>det.outer_rad;
     else if ((int)line.find("outer_depth")!=-1) sin>>det.outer_depth;

     

  }

  
}
// ########################################################
// CC neutrino cross-section (cm2) - various models fitted
// ########################################################

double dsigGR(double E, int type, int AntiNu)
{
  //check if it is anti nuE
  if(type!=12 || AntiNu!=-1) return 0.;

  double GF = 1.166E-5; //GeV^-2
  
  double A=3.4986E-27;
  double gammaW=2.085; //GeV/c^2

  double cross_section_num=A*4*GF*GF*E*mW*mW*mW*mW*me;
  double cross_section_denom=3*2*PI*((mW*mW-2*me*E)*(mW*mW-2*me*E)+mW*mW*gammaW*gammaW);
  if(cross_section_denom==0){cout<<"GR cross section /0 error"<<endl;}
  double cs=cross_section_num/cross_section_denom;

  return cs;
}
double dsigCC(double E, int CCmode, int type,int AntiNu )
{
  
  double f=0.;
  double p[4];
  
  AntiNu=int(abs(AntiNu-1)/2);//return AntiNu to be like a bool
  // The value below determines when we switch from the parameterizations 
  // of the neutrino cross sections at ultra-high energis (e.g. CTTW standard values) 
  // to the cross sections to the Ghandi parameterization. These transitions were determined
  // using the parameterization made for this code. They are likely a bit different if the user 
  // switches the cross section to the upper or lower cross section models. 

  double E_switch = 2.00e6;  // GeV

  // If the energy is below E_sigma_switch, set CCmode to the Ghandi cross-section.
  // If the particle is a neutrino, AntiNu = 0 and the cross-section is set to the
  // the Ghandi model for neutrinos. If it is an anti-neutrino, AntiNu=1 and the
  // cross-section is set to the Ghdni model for anti-neutrinos.

  if( E < E_switch ){
  CCmode = 3 + AntiNu;
  } 

  // Connolly+, 2011 middle model (ARW's parametrization)
  double p0[4] = { -5.35400180e+01,   2.65901551e+00, -1.14017685e-01,   1.82495442e-03};
  // Connolly+, 2011 lower model (ARW's parametrization)
  double p1[4] = {-4.26355014e+01,   4.89151126e-01,   2.94975025e-02,  -1.32969832e-03};
  // Connolly+, 2011 upper model (ARW's parametrization)
  double p2[4] = {-5.31078363e+01,   2.72995742e+00,  -1.28808188e-01,   2.36800261e-03};

  // Gandhi, Quigg, Reno 1995 Neutrino cross section
  double p3[4] = { -6.24043607e+01,   4.21769574e+00, -2.06814586e-01,   3.70730061e-03};
  // Gandhi, Quigg, Reno 1995 Anti-Neutrino cross section
  double p4[4] = { -6.43574494e+01,   4.41740442e+00, -2.10856220-01,   3.65724741e-03};

  double log10_E_eV = log10(E)+9.;
  for (int ii = 0 ; ii<4; ii++)
  {
    if(CCmode==0) p[ii] = p0[ii];
    if(CCmode==1) p[ii] = p1[ii];
    if(CCmode==2) p[ii] = p2[ii];
    if(CCmode==3) p[ii] = p3[ii];
    if(CCmode==4) p[ii] = p4[ii];

    f += p[ii]*pow(log10_E_eV,ii);
  }

  f = pow(10,f);
  return f;

    
        
  // 	double l1=log10(E);
  // 	double l2=l1*l1;
  // 	double l3=l2*l1;
  // 	double l4=l3*l1;
  // 	double l5=l4*l1;
  // 	double l6=l5*l1;
  // 	double l7=l6*l1;
  // 	f=pCC0+pCC1*l1+pCC2*l2+pCC3*l3+pCC4*l4+pCC5*l5+pCC6*l6+pCC7*l7;
  
  //      f = 6.37994*pow(E,0.355991)*1e-36;
  
  /* CKMT */
  //      f = (-36.3345965603+7.14693605311*pow(E,0.293313250614))*1.e-36;
  
  /* ALLM */
  // H. Abramowicz et al., Phys. Lett. B 269, 465 (1991);
  // H. Abramowicz and A. Levy, hep-ph/9712415.
  //	f = (-280.544665122+10.3452620208*pow(E,0.317119535055))*1.e-36;
  
  /* ASW */   // Saturation of pdfs
              // N. Armesto et al., Phys. Rev. D 77, 013001 (2008).
              // N. Armesto et al., Phys. Rev. Lett. 94, 022002 (2005).
              //      f = (-799.252409182+52.4932827684*pow(E,0.244551044541))*1.e-36;
  
  /* Sarkar */  // Default model used in Auger
                  // A. Cooper-Sarkar and S. Sarkar, JHEP 0801, 075 (2008).
                  // Amanda Cooper-Sarkar, Philipp Mertsch, Subir Sarkar. JHEP 08, 042 (2011).
                  //      f = (-649.265343982+26.4437052803*pow(E,0.296160447336))*1.e-36;
  
  // Sarkar model (Yann's parametrization)
  //    double AS=-0.391641;
  //    double BS=0.635232;
  //    double CS=-0.0158144;
  //    f= (pow(10,AS+BS*log10(E)+CS*pow(log10(E),2)))*1.e-36;
        
      
}


double dsigNC(double E, int CCmode, int type,int AntiNu)
{
  double f=0.; 
  
      
  AntiNu=int(abs(AntiNu-1)/2);//return anti ness as a bool
  double p[4];

  // The value below determines when we switch from the parameterizations 
  // of the neutrino cross sections at ultra-high energis (e.g. CTTW standard values) 
  // to the cross sections to the Ghandi parameterization. These transitions were determined
  // using the parameterization made for this code. They are likely a bit different if the user 
  // switches the cross section to the upper or lower cross section models. 
  
  double E_switch = 2.00e6;  // GeV

  // If the energy is below E_sigma_switch, set CCmode to the Ghandi cross-section.
  // If the particle is a neutrino, AntiNu = 0 and the cross-section is set to the
  // the Ghandi model for neutrinos. If it is an anti-neutrino, AntiNu=1 and the
  // cross-section is set to the Ghdni model for anti-neutrinos.

  if( E < E_switch ){
  CCmode = 3 + AntiNu;
  }

  // Connolly+, 2011 middle model (ARW's parametrization)
  double p0[4] = { -5.41463399e+01,   2.65465169e+00,  -1.11848922e-01,   1.75469643e-03};
  // Connolly+, 2011 lower model (ARW's parametrization)
  double p1[4] = {-4.42377028e+01, 7.07758518e-01, 1.55925146e-02, -1.02484763e-03};
  // Connolly+, 2011 upper model (ARW's parametrization)
  double p2[4] = {-5.36713302e+01,   2.72528813e+00,  -1.27067769e-01,   2.31235293e-03};
      
  // Gandhi, Quigg, Reno 1995 Neutrino cross section
  double p3[4] = { -6.33753554e+01,   4.26790713e+00,  -2.07426844e-01,   3.68501726e-03};
  // Gandhi, Quigg, Reno 1995 Anti-Neutrino cross section
  double p4[4] = { -6.33697437e+01,   4.11592385e+00,  -1.90600183e-01,   3.22478095e-03};

  double log10_E_eV = log10(E)+9.;
  for (int ii = 0 ; ii<4; ii++){
    if(CCmode==0) p[ii] = p0[ii];
    if(CCmode==1) p[ii] = p1[ii];
    if(CCmode==2) p[ii] = p2[ii];
    if(CCmode==3) p[ii] = p3[ii];
    if(CCmode==4) p[ii] = p4[ii];

    f += p[ii]*pow(log10_E_eV,ii);
  }

  f = pow(10,f);
  
  return f;
}

// ###################################################
// 1./(gamma*c*muon0) with muon0 lifetime of muon (cm^-1)
// ###################################################
double dPdesdx(double E, int type)
{
    double f;
    if(type==13||type==14) f=mmuon/(E*muondl);
    if(type==15||type==16) f=mtau/(E*taudl);
        
    return f;
}

double elost(double E, double dens, int ELOSSmode,int type)
{
  double f;
  //double z = 0.;
  //  dE/dX =      E*beta(E)      +     alpha(E)
  
  // this is a super-kludgy way to account for beta being different for iron <A>=56.84, <Z>=26; rock <A>=22, <Z>=11; and water <A>=11.9, <Z>=6.6
  // material properties are not tracked in this simulation but densities are.

  int lyr = 0; // initialize to iron
  if (dens < 7.75) lyr = 1; // density jump between outer core and mantle
  if (dens < 2.0) lyr = 2;  // density jump between rock and water 
  if(type==13)
  {
  double factor[3] = {0.9304, 1.0, 1.1092}; // ratio Z/A for iron, rock, and water divided by Z/A=0.5 for rock   

  // The correction below was based on the claim that the photonuclear energy loss is proportional to <A> as well as the density in Palomares-Ruiz, Irimia, & Weiler, Phys. Rev. D, 73, 083003 (2006)
  // Searching through the references, this claim is demonstrably false. See S. I. Dutta, M. H. Reno, I. Sarcevic, and D. Seckel, Phys. Rev. D 63, 094020 (2001).
  // Earlier runs of the code used the line below but it has been commented out.
  // if(dens<1.1) factor = 0.55; // This kluge corrects for the change of <A>=22 in rocks vs <A>=12 of H2O
  //printf("factor %1.2f \n", factor);
  
  //f = E * beta9fit(&E,&lyr,ELOSSmode) + mfuncalph(&E, &lyr,type);
  f = 2e-3*factor[lyr] + E * beta9fit(&E,&lyr,ELOSSmode,type);
	
	//f = E*(emlost->Eval(E)) + funca->Eval(E);
  // cout << "\tE " << E << endl;
  //cout << "\tlyr " << lyr << " dens = " << dens << endl;
  // cout << " f = " << f << endl;
  // cout << " E * beta9fit(&E,&lyr) + mfuncalph(&E, &lyr); " << E * beta9fit(&E,&lyr) + mfuncalph(&E, &lyr) << endl << endl;
  // cout << " mfuncalph(&E,&lyr) " << E << " " << lyr << " " << mfuncalph(&E, &zlyr << endl << endl;
  //cout << " beta9fit(&E,&lyr) " << E << " " << beta9fit(&E,&lyr,0) << " " << beta9fit(&E,&lyr,1) << endl << endl;
  }
  if(type==15)
  {
    f = E * beta9fit(&E,&lyr,ELOSSmode,type) + funcalph(&E, &lyr,type);
    //f = E*(emlost->Eval(E)) + funca->Eval(E);
    // cout << "\tE " << E << endl;
    //cout << "\tlyr " << lyr << " dens = " << dens << endl;
    // cout << " f = " << f << endl;
    // cout << " E * beta9fit(&E,&lyr) + funcalph(&E, &lyr); " << E * beta9fit(&E,&lyr) + funcalph(&E, &lyr) << endl << endl;
    // cout << " funcalph(&E,&lyr) " << E << " " << lyr << " " << funcalph(&E, &zlyr << endl << endl;
    //cout << " beta9fit(&E,&lyr) " << E << " " << beta9fit(&E,&lyr,0) << " " << beta9fit(&E,&lyr,1) << endl << endl;
  }
  if(f==0){cout<<"elost is zero"<<endl;}
  return f;
}

// ###################################################
// ###################################################
double funcalph(double *x, int *par, int type)
// double tfuncalph(double *x)
{ 
  double f;
  if(type==15)
  {
    
    double p=sqrt(x[0]*x[0]-mtau2);
    double b=p/x[0];
    double b2=b*b;
    double gamma=x[0]/mtau;
    double EE=Cbb2*p*p/(me2+mtau2+Cbb2*x[0]);
    double X=log10(b*gamma);
    double factor[3] = {0.9304, 1.0, 1.1092}; // ratio Z/A for iron, rock, and water divided by Z/A=0.5 for rock 
    
    f=Cbb1/(b2)*(log(Cbb2*b2*gamma*gamma/I2)-2*b2+EE*EE/(4*x[0]*x[0])-delta(X));
    f *= factor[par[0]];
    //cout << "\tlyr " << par[0] << " factor " << factor[par[0]] << " val " << f << endl; 
  
  }
  if(type==13)
  {
    double f;
    double p=sqrt(x[0]*x[0]-mmuon2);
    double b=p/x[0];
    double b2=b*b;
    double gamma=x[0]/mmuon;
    double EE=Cbb2*p*p/(me2+mmuon2+Cbb2*x[0]);
    double X=log10(b*gamma);
    double factor[3] = {0.9304, 1.0, 1.1092}; // ratio Z/A for iron, rock, and water divided by Z/A=0.5 for rock 
    
    f=Cbb1/(b2)*(log(Cbb2*b2*gamma*gamma/I2)-2*b2+EE*EE/(4*x[0]*x[0])-delta(X));
    f *= factor[par[0]];
    //cout << "\tlyr " << par[0] << " factor " << factor[par[0]] << " val " << f << endl; 
    
  }
  return f;
  
}

// ###################################################
// ###################################################
double delta(double X)
{
  double f=0;
  
  if(X > X1)
  {
    f=4.6052*X+CC;
  }
  else if((X>X0)&&(X<X1))
  {
    f=4.6052*X+CC+aa*pow((X1-X),mm);
  }
  
  return f;
}

// ###################################################
// beta(E) function in dE/dX - various models fitted
// ###################################################
//mac double beta9fit(double *x, double *par)
//double beta9fit(double *x)
double beta9fit(double *x, int *par, int ELOSSmode, int type)
{
    if(type==15||type==16) //tau
    {
        //double f=0.;
        double b0 = 0.;
        double b1 = 0.;
        double b2 = 0.;
        
        /* ALLM */
        if(ELOSSmode==0)
        {
        //b0 = 2.05820774222e-07;
        //b1 = 4.93367455295e-09;
        //b2 = 0.227781737887;
        b0 = -7.78527765e+00;
        b1 = -2.80672147e-02;  
        b2 = 6.38891661e-03;
        //printf("ALLM \n");
        }
        /* CKMT */
        /*
        double bo = 1.93693562238e-07;
        double b1 =  4.6247859382e-09;
        double b2 =  0.224891839933;
        */
        
        /* ALLM */
        // H. Abramowicz et al., Phys. Lett. B 269, 465 (1991);
        // H. Abramowicz and A. Levy, hep-ph/9712415.
        /*
        double b0=2.05820774222e-07;
        double b1=4.93367455295e-09;
        double b2=0.227781737887;
        */
        
        /* ASW */  // Saturation of pdfs
                    // N. Armesto et al., Phys. Rev. D 77, 013001 (2008).
                    // N. Armesto et al., Phys. Rev. Lett. 94, 022002 (2005).
        if(ELOSSmode==1)
        {
        b0=-4.77043758142e-08;
        b1=1.9031520827e-07;
        b2=0.0469916563971;
        b0 = -8.61998129e+00;
        b1 =  1.57820040e-01; 
        b2 = -2.24340096e-03;
        //printf("ASW \n");
        }
        double log10E = log10(x[0]);
        double f_phot = pow(10., b0 + b1*(log10E+9.) + b2*(log10E+9.)*(log10E+9.));
        //f=b0+b1*pow(x[0],b2);
        //printf("%1.2e \n", f);
        
        double f_brem = pBrem[par[0]][0]*pow(1.-exp(-pow((log10E)/pBrem[par[0]][1], pBrem[par[0]][2])),pBrem[par[0]][3]);
        double f_pair = pPair[par[0]][0]*pow(1.-exp(-pow((log10E)/pPair[par[0]][1], pPair[par[0]][2])),pPair[par[0]][3]);

        // cout << "\tpar[0] " << par[0] << endl;
        // cout << "\tBrem " << pBrem[par[0]][0] << " " <<  pBrem[par[0]][1] << " " <<  pBrem[par[0]][2] << " " <<  pBrem[par[0]][3] << endl;
        // cout << "\tPair " << pPair[par[0]][0] << " " <<  pPair[par[0]][1] << " " <<  pPair[par[0]][2] << " " <<  pPair[par[0]][3] << endl;
        // cout << "\ttest " << 1.-exp(pow(log10(x[0])/pBrem[par[0]][1], pBrem[par[0]][2])) << endl;
        // cout << "\tpair, brem, photoN " << f_pair << " " << f_brem << " " << f_phot << endl;
        return f_pair + f_brem + f_phot;

    };
    if(type==13||type==14) //muon 
    {

            /*Energy loss parameters here are the sum ofbremmstrahlung, pair production, and photonuclear interactions*/
            /*Photonuclear losses are characterized using the 2 models below*/
            double f=0.;
            double b0 = 0.;
            double b1 = 0.;
            double b2 = 0.;
            double b3 = 0.;

            int lyr = par[0];
                
            /* BB */
            /* Bezrukov and Bugaev Model for photonuclear losses*/
            /* L. B. Bezrukov and E. V. Bugaev, Sov. J. Nucl. Phys. 33, 635 (1981). */
            if(ELOSSmode==1)
            {
                if (lyr==1){//Rock
                b0 = 5.83467127e-7;
                b1 = 1.43408350e-6;  
                b2 = -1.68421469e-7;
                b3 = 7.00793338e-9;
                //printf("Rock \n");
                }  
            if(lyr==2){//Water
                b0 = 5.12473265e-7;
                b1 = 9.99525050e-7;
                b2 = -1.11322743e-7;
                b3 = 4.72858936e-9;
                //printf("Water \n");
                }    
            if(lyr==0){//Iron
                b0 = 8.74681421e-7;
                b1 = 2.93621999e-6;
                b2 = -3.62818552e-7;
                b3 = 1.47773834e-8;
                //printf("Iron \n");
                }
                //printf("BB \n");
            }

            /* ALLM */
            /* ALLM Model for photonuclear losses*/
            /* S. Dutta, M. H. Reno, I. Sarcevic, D. Seckel Phys.Rev. D63 (2001) 094020 */
            if(ELOSSmode==0)
            {
                if (lyr==1){//Rock
                b0 = 8.75384980e-8;
                b1 = 2.05815984e-6;  
                b2 = -3.38061300e-7;
                b3 = 2.00850609e-8;
                //printf("Rock \n");
                }
                if(lyr==2){//Water
                b0 = -3.08528426e-8;
                b1 = 1.67990356e-6;
                b2 = -2.94785318e-7;
                b3 = 1.87902903e-8;
                //printf("Water \n");
                }
                if(lyr==0){//Iron
                b0 = 4.49017534e-7;
                b1 = 3.48261510e-6;
                b2 = -5.13196702e-7;
                b3 = 2.64637027e-8;
                //printf("Iron \n");
                }
                //printf("ALLM \n");
            }
            double log10E = log10(x[0]);
            f = b0+b1*log10E+b2*log10E*log10E+b3*log10E*log10E*log10E;

            return f;
    }
    
  
}



// #######################################################
// Mean Earth density along a chord at angle theta
// #######################################################
double mean_dens_chord(double theta)
{
  double f;
  double z= 0.;
  f = mean(&theta, &z);
  // f = meandens->Eval(theta);
  // cout << "f = meandens->Eval(theta);" << meandens->Eval(theta) << endl;
  // cout << "f = mean(theta);" << mean(&theta, &z) << endl;
  
  return f;
}


//// ###################################################
//// Mean Earth density along a chord of length Lmax
//// ###################################################
////mac double mean(double *x, double *par)
double mean(double *x, double *par)
{
  // double f;
  
  double Lmax=2*R0*cos(PI-PI*x[0]/180.);
  
  //TF1 *densi = new TF1("dens",earthdens,0.,Lmax,1); // (name, function, xmin, xmax, number of parms.
  //densi->SetParameters(Lmax,0); // set parameters to Lmax (par[0] = Lmax)
  
  //f = (densi->Integral(0.,Lmax))/Lmax;
  //cout << " f = " << f << endl;
  
  // stupid numerical integral
  double sum = 0.;
  int N = 1000000;
  double dx = Lmax/((double) (N));
  for (int ii=0; ii<=N; ii++)
  {
  double x_val = ((double) ii)*dx;
  sum += earthdens(&x_val, &Lmax) * dx;
  }
  sum /= Lmax;
  //cout << " sum = " << sum << endl;
  
  //delete densi;
  
  return sum;
}

// #######################################################
// Earth density as a function of radius - read from table
// #######################################################
double earthdens( double *x, double *par)
{
  double f;
  
  double Current_Radius = sqrt(R02 - (x[0]*par[0])+x[0]*x[0]);
  
  f = terra->GetDensity(Current_Radius);
  //printf("check: %1.2e %1.2e\n", Current_Radius, f);
  
  return f;
}				

//initializes the reaction data arrays from the pythia tables
void initialize_reaction(int tau_type[][6], int mu_type[][6],double tau_ene[][6],double mu_ene[][6])
{
  
    string line;
    string list_of[100000];
    string reaction;
    string particle[6];
    string energy[6];
    string delimiter_p=":";
    string delimiter_e=",";
    string products[6];
    int counter;
    int done;

    //taus
    ifstream tau_decays("pythia_tau.txt");
    
    for(int count=0;count<100000;count++)
    {
        getline(tau_decays,line);
        list_of[count]=line;
       
    }
    
    //parse string for products and fractional energies.
    for(int i=0;i<100000;i++)
    {
        for(int z=0;z<6;z++) {products[z]="0,0.0";};
        counter=0;
        done=0;
        
        for(int j=0;j<6;j++)
        {
            if(list_of[i].find(delimiter_p)!=string::npos)
            {
                products[j]=list_of[i].substr(0,list_of[i].find(delimiter_p));
                list_of[i].erase(0,list_of[i].find(delimiter_p)+1);
                counter++;
            }
            if(list_of[i].find(delimiter_p)==string::npos && done==0)
            {
                products[counter]=list_of[i];
                done=1;
            }
        }
        
        for(int s=0;s<6;s++)
        {
            tau_type[i][s]=convert_types(stoi(products[s].substr(0,products[s].find(delimiter_e))));
            products[s].erase(0,products[s].find(delimiter_e)+1);
            tau_ene[i][s]=stod(products[s]);
        }
    }
    tau_decays.close();

    //muons
    ifstream mu_decays("pythia_muon.txt");
    
    for(int count=0;count<10000;count++)
    {
        getline(mu_decays,line);
        list_of[count]=line;
    }

    //parse string for products and fractional energies.
    for(int i=0;i<100000;i++)
    {
        for(int z=0;z<6;z++) {products[z]="0,0.0";}
        counter=0;
        done=0;
       for(int j=0;j<6;j++)
        {
            if(list_of[i].find(delimiter_p)!=string::npos)
            {
                products[j]=list_of[i].substr(0,list_of[i].find(delimiter_p));
                list_of[i].erase(0,list_of[i].find(delimiter_p)+1);
                counter++;
            }
            if(list_of[i].find(delimiter_p)==string::npos && done==0)
            {
                products[counter]=list_of[i];
                done=1;
            }
        }
        
        for(int s=0;s<6;s++)
        {
            mu_type[i][s]=convert_types(stoi(products[s].substr(0,products[s].find(delimiter_e))));
            products[s].erase(0,products[s].find(delimiter_e)+1);
            mu_ene[i][s]=stod(products[s]);
        }
    }
    mu_decays.close();
}
void make_dirs(string dirs)
{
  mkdir(dirs.c_str(),0755);
  mkdir((dirs+"/particles").c_str(),0755);
  mkdir((dirs+"/events").c_str(),0755);
  mkdir((dirs+"/LUT").c_str(),0755);
}

int convert_types(int pythia_type) 
{
    int type=0;
      switch (pythia_type){

        case -15:
            type=-15;
            break;
        case 15:
            type=15;
            break;
        case -16:
            type=-16;
            break;
        case 16:
            type=16;
            break;
        case 13:
            type=13;
            break;
        case -13:
            type=-13;
            break;
        case 14:
            type=14;
            break;
        case -14:
            type=-14;
            break;
        case 12:
            type=12;
            break;
        case -12:
            type=-12;
            break;
        case 11:
            type=11;
            break;
        case -11:
            type=-11;
            break;    
        default:
            type=0;   
    }

    /*
    switch (pythia_type){

        case -15:
            type=-6;
            break;
        case 15:
            type=6;
            break;
        case -16:
            type=-3;
            break;
        case 16:
            type=3;
            break;
        case 13:
            type=5;
            break;
        case -13:
            type=-5;
            break;
        case 14:
            type=2;
            break;
        case -14:
            type=-2;
            break;
        case 12:
            type=1;
            break;
        case -12:
            type=-1;
            break;
        case 11:
            type=4;
            break;
        case -11:
            type=-4;
            break;    
        default:
            type=0;   
    }
    */
    return type;
}



