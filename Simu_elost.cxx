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


//updated on 1/8/24 by Ryan Krebs

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

using namespace std;
/*
class prop_particle
{
  public:  
    int ID;
    double energy;
    double position,x,y,z,u,v,w;
    double type;
    int NC,DC,CC,GR;

    void set_pos(double temp_x, double temp_y, double temp_z);
    void set_traj(double temp_u, double temp_v, double temp_w);
    void step_pos(double length);
    void reset_class();

};  
void prop_particle::set_pos(double temp_x, double temp_y, double temp_z)
{
  x=temp_x;
  y=temp_y;
  z=temp_z;
}
void prop_particle::set_traj(double temp_u, double temp_v, double temp_w)
{
  u=temp_u;
  v=temp_v;
  w=temp_w;
}
void prop_particle::step_pos(double length)
{
  x=length*u+x;
  y=length*v+y;
  z=length*w+z;
}
*/
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
/*
class charged_lep
{
  public:
    double xi,xf,yi,yf,zi,zf,Ei,Ef;
    int p_type,anti;
  

    void fill_initial(double x,double y,double z, double E,int antiness, int p);
    void fill_end(double x,double y,double z, double E);
    void empty_class();
    void save_lep(ofstream *output_file,int code,int t_num,int n_num);
    void print_values();
};
void charged_lep::print_values()
{
  printf("p_type: %i, anti: %i,\nxi: %f, yi: %f, zi: %f, Ei: %f,\nxf: %f, yf: %f, zf:%f, Ef: %f \n "
  ,p_type,anti,xi,yi,zi,Ei,xf,yf,zf,Ef);
}
void charged_lep::fill_initial(double x,double y,double z, double E, int antiness,int p)
{
  xi=x;
  yi=y;
  zi=z;
  Ei=E;
  p_type=p;
  anti=antiness;
}
void charged_lep::fill_end(double x,double y,double z, double E)
{
  xf=x;
  yf=y;
  zf=z;
  Ef=E;
}
void charged_lep::empty_class()
{
  //printf("emptying charged lepton\n");
  xi=0;
  yi=0;
  zi=0;
  xf=0;
  yf=0;
  zf=0;
  Ei=0;
  Ef=0;
  p_type=0;
}
void charged_lep::save_lep(ofstream *output_file,int code,int t_num,int n_num)
{
  if(Ei!=0&&Ef!=0)
  {
    *output_file<<xi<<","<<yi<<","<<zi
  <<","<<xf<<","<<yf<<","<<zf
  <<","<<Ei<<","<<Ef<<","<<p_type*anti
  <<","<<code<<","<<t_num<<","<<n_num<<endl;
  
  //printf("saving");
  }

}
*/
typedef struct {
  stack<int> parent_part_ind;
  stack<int> part_type;         // pythia code - only positive values
  stack<double> part_energy;    // E_min<part_energy<E_init
  stack<double> part_pos;       // position inside earth
  stack<int> anti;              // true if anti particle
  stack<int> generation;        // number of interactions before this particle was produced
  stack<int> NC_num;            // number of NC preceding creation of this particle
  stack<int> CC_num;            // number of CC preceding creation of this particle
  stack<int> GR_num;            // number of GR preceding creation of this particle
  stack<int> dc_num;            // number of decays it took to get to the current particle
  stack<double> x_pos;          // x coord. where part. is created
  stack<double> y_pos;          // y coord. where part. is created
  stack<double> z_pos;          // z coord. where part. is created
  stack<double> traversed_gram; // traversed grammage up until particle is created 
  stack<int> start_in_volume;

} particle_info_def;
typedef struct 
{
  double xi,yi,zi;
  double xf,yf,zf;
  double Ei,Ef;
  int p_type;

} charged_leptons;


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
  bool save_sto;
  int n_throws;
  int n_traj;
  bool save_nu_events;
  bool save_sto_events;
  bool save_dec_events;
  bool save_emerging;
  bool use_sto_inst_of_cont;
}config_init;   // data struct to hold the value read from config file



typedef struct
{
  //last 6 are on the surface. first 6 are on a set interface --not used
  double xi;
  double yi;
  double zi;
  double xf;
  double yf;
  double zf;
  double xe;    // x coord of initial earth point
  double ye;    // y coord of initial earth point
  double ze;    // z coord of initial earth point
  double eex;   // x coord of exit earth point
  double eey;   // y coord of exit earth point
  double eez;   // z coord of exit earth point
  //int input_types;
  //double input_e;
}input_file;//input file parameters setting startpoint, end point, and etrance point to the volume


typedef struct 
{
  //cylindircal
  double inner_rad;
  double outer_rad;
  double inner_depth;
  double outer_depth;
  //spherical
  double inner_sphere;
  double outer_sphere;

  //gen
  double traj_weights;

}det_geom;


//MYEVT_DEF event;
particle_info_def particle_data;
reaction_tables_def reaction_data;
config_init config;
det_geom det;
double decay_length(double P, double E, int particle_type);
bool in_volume(double x,double y,double z);
//void generate_events(det* detec);

string make_particle_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp);
string make_event_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp, int p_type);
string make_lepton_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp, int p_type);
// initialize reaction from pythia table and convert pythia type to tag code used in code
void initialize_reaction(int tau_type[][6], int mu_type[][6],double tau_ene[][6],double mu_ene[][6]);

// converts pythia tags to tags in this code. converts anti particles to normal matter
int convert_types(int pythia_type);

//load values from config file
void load_config();
void load_geo();
int load_input(input_file *in,int length, string filename);

void set_points_from_angle(input_file *input,double angle);
// -------------------------------------------------
// For lepton energy loss: dE/dX = -alpha + beta(E)*E 
//double funcalph(double *x, int *par, int type);

double delta(double X);

// Parameterisations for beta 
//double beta9fit(double *x, int *par, int ELOSSmode, int type);	

// Elost by tau and muon dE/dX in GeV/(g/cm^2)
//double elost(double E, double dens, int ELOSSmode, int type);

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
  int taus_passed_through=0;
  int decay_num=0;
  //charged_lep lep;
  continuous_loss_prop cont;
  stochastic_lepton_prop sto;
  sto.load_tables();
  // call function to load config
  load_config(); 
  make_dirs(config.data_dir);
  load_geo();//might not be necessary now that points are defgined outside of sim
  //charged_leptons temp_lep;
  double tau_x = 1;
  double tau_y=1;
  double tau_z=1;
  string in_file="in_files/59068000000_Example_Trajectories.csv";
  int input_num = config.n_traj;
  if(config.detector==0)input_num=1;
  double traj_weight=1/59068000000;
  input_file *input = new input_file[input_num];
  //new input_file input[input_num];
  
  //load_input(input,input_num,in_file);
  
  
  // initializes arrays to hold decay products and populates them from pythia file
  for(int i=0;i<100000;i++){for(int j=0;j<6;j++){reaction_data.tau_type[i][j]=reaction_data.mu_type[i][j]=0;reaction_data.tau_energy[i][j]=reaction_data.mu_energy[i][j]=0.0;};}  
  initialize_reaction(reaction_data.tau_type,reaction_data.mu_type,reaction_data.tau_energy,reaction_data.mu_energy);
//{reaction_data.tau_type[i][j]=0;reaction_data.mu_type[i][j]=0;reaction_data.tau_energy[i][j]=0.0;reaction_data.mu_energy[i][j]=0.0;};}  
  int gr_counter=0;
  //cout<<config.detector<<" "<<config.det_volume<<endl;
  

  int type_to_save[5]={-1,15,-1,-1,-1}; //max of 5 types since electrons are forgotten
 //hard coding this bc I dont feel like adding it as an option
  //if(config.save_charged){type_to_save[0]=13;type_to_save[1]=15;} //add muons and taus to pareticle type to savce
  //if(config.save_neutrinos){type_to_save[2]=12;type_to_save[3]=14;type_to_save[4]=16;} //add neutrinos to type to save
  
  //for(int i=0;i<5;i++) cout<<type_to_save[i]<<" "; //print which particles will be saved
  //cout<<endl;

  cout << "Lepton Propagation code" << endl;

  
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

  if(config.detector==0)set_points_from_angle(input,angle);
  if(config.detector==1)load_input(input,input_num,in_file);
  
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
    
  //cout<<config.starting_type<<","<<atoi(argv[8])<<endl;
  //if(argc>8) config.starting_type=atoi(argv[8]);
  cout<<config.starting_type<<endl;
  //cout<<angs_temp<<endl;
  //-------------------------------------------------
  // Output file names using input arguments
  string nameEnergies="";
  string nameEvents="";
  string nameOutLep="";
  nameEnergies=make_particle_dir(argc,argv,config.data_dir,es_temp,angs_temp); 
  nameEvents=make_event_dir(argc,argv,config.data_dir,es_temp,angs_temp,config.starting_type); 
  //nameOutLep=make_lepton_dir(argc,argv,config.data_dir,es_temp,angs_temp,config.starting_type);
  


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
  outEnergies << "type, anti, NC, CC, GR, DC, Gen, InitNuNum, InitNeutrinoType, OutEnergy, InitEnergy, Part_Pos.\n";
  //cout<<nameEnergies<<endl;
  ofstream outEvents(nameEvents.c_str());
  outEvents<<"vert_x,vert_y,vert_z,tra_x,tra_y,tra_z,E_nu,inel,part_type,i_type,had_or_em,nc_num,dc_num,gr_num,traj_num,p_thrown,sto_index"<< setprecision(9)<<endl;

  

  //ofstream outLep(nameOutLep.c_str());
  //outLep<< setprecision(9)<<"logEi,f,logEd,x,y,z,vx,vy,vz,sto_type,traj,p_num,p_type,sto_num"<<endl;

  //string outGrammage="lepton_test/grammage_"+es_temp+"_"+to_string(config.starting_type)+".dat";

  //ofstream out_gram(outGrammage.c_str());
  //out_gram<<"t_num,n_num,km,grammage,event"<<endl;
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
  //cout << "Average Earth density " << mean_dens_chord(refTheta) << " g/cm^3" << endl;
  //cout << endl;
  
  cout << "======================================" << endl;
  //cout << "Emerging Lepton Energy file: " << endl;
  //cout << nameEnergies << endl;
  //cout << endl;
  
  //cout << "======================================" << endl;
  

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


  for( int i=0;i<input_num;i++)//loop over trajectories
  { 
    if(i%10000==0)cout<<"Did "<<i<<" trajectories..."<<endl;
    //cout<<i<<endl;
    int num_count=0;
    //---------------------------------
    double sum_grammage = 0.;
    
    double check_grammage=0;
    //cout<<input[i].eex<<","<<input[i].xi<<"..."<<input[i].eez<<","<<input[i].zi<<endl;
    double maxL=sqrt((input[i].eex-input[i].xi)*(input[i].eex-input[i].xi)
                    +(input[i].eey-input[i].yi)*(input[i].eey-input[i].yi)
                    +(input[i].eez-input[i].zi)*(input[i].eez-input[i].zi));
    //cout<<maxL<<endl;
    if(maxL<10)
    {
      cout<<"didnt load correctly";
      return 0;
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
    temp_pos[0]=temp_pos[0]+dl*x_step;
    temp_pos[1]=temp_pos[1]+dl*y_step;
    temp_pos[2]=temp_pos[2]+dl*z_step;
    //printf("*** ii %d %1.5f %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
    //if(ii%100000 ==0) printf("ii %d %1.2e %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
    }


    //save charged leptons entering the volume!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    //----------------------------------
    //cout<<"primary trajectory: "<<i<<endl;
    bool got_event=false;//!got_event && num_count<1000 line below
    int out_leptons=0;

    //while(num_count<config.n_throws)//!got_event for running until getting 1 event or num_count<N for throwing N particles
    while(out_leptons<500)
    {
    //cout<<num_count<<endl;
    //if(num_count==10000)cout<<"10,000...";
    //if(num_count==100000)cout<<"100,000..."<<endl;
    //if(num_count%1000000==0)cout<<"did "<<num_count<<" thrown particles"<<endl;
    //if(out_leptons%100==0)cout<<"caught "<<out_leptons<<" taus"<<endl;
    //cout<<num_count<<endl;
    //if(config.save_sto)lep.empty_class();
    num_count++;
    //temp_lep.Ei=0;
    //temp_lep.Ef=0;
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
    int anti=1;
    anti=config.anti;
    if(config.starting_type==16) {part_type=16; anti=1;}
    if(config.starting_type==14) {part_type=14; anti=1;}
    if(config.starting_type==12) {part_type=12; anti=1;}
    if(config.starting_type==-12) {part_type=12; anti=-1;}
    if(config.starting_type==-14) {part_type=14; anti=-1;}
    if(config.starting_type==-16) {part_type=16; anti=-1;}

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
    
    int generation=0;
    int NC_num=0;
    int CC_num=0;
    int GR_num=0;
    int dc_num=0;
    

    //stuff for w decays

    const double e_rest=(mW*mW+me2)/(2*mW);
    const double m_rest=(mW*mW+mmuon2)/(2*mW);
    const double t_rest=(mW*mW+mtau2)/(2*mW);
    
  
    traversed_grammage = 0.; //initialize traversed grammage for each event
    double traversed_path=0.;

    
    
    //======================= Loop over stacks until empty
    do 
    { 
      
      has_been_saved=false;
      bool save_cond=false;
      bool exit_cond=false;
      int start_in_volume=0;
      bool save_lep=false;
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
        if(config.save_sto)start_in_volume=particle_data.start_in_volume.top();
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
        if(config.save_sto)particle_data.start_in_volume.pop();
      }
      double tau_path=0;
      //if(loop_num!=0)cout<<"stack loop of type "<<part_type<<" and energy "<<part_energy<<"at x,y "<<pos[0]<<","<<pos[1]<<endl;
      if(part_type==11) continue; //ignores particles of electron flavor
      if(part_energy<Elim) continue; //ignore particles below threshold in case they make it through
      if(config.save_sto&&start_in_volume==1 && (part_type==13 || part_type==15))
      {
        //empty_class();
        //lep.fill_initial(pos[0],pos[1],pos[2],part_energy,anti,part_type);
        //temp_lep.xi=pos[0];
        //temp_lep.yi=pos[1];
        //temp_lep.zi=pos[2];
        //temp_lep.Ei=part_energy;
        //temp_lep.p_type=part_type;
        //temp_lep.Ef=0;
      }
      bool checked_in=false;
      //cout<<part_type<<","<<anti<<endl;

      // Flag to see how fast code is running
      //if(!((float)i/100000-(int)(i/100000))) cout<< i << endl;
      

      //brkcnt=0;
      //prop_mode =0;
      
      //======================= Start propagation along chord of length Lmax
      //printf("Ldist, Lmax %1.2e %1.2e\n", Ldist, Lmax);
      //Ldist=0.;
      int sto_num=0;      
      bool broken=false;
      bool left_volume=false;
      bool entered_volume=false;
      int sto_index=0;
      //while(part_pos<maxL &&  !left_volume) 
      while(part_pos<maxL)
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
        
        //cout<<"dist-";
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
            pos[0]=pos[0]+(after_step-before_step)*x_step;
            pos[1]=pos[1]+(after_step-before_step)*y_step;
            pos[2]=pos[2]+(after_step-before_step)*z_step;
          }

          
          if(in_volume(pos[0],pos[1],pos[2])&& !in_vol)
          {
            entered_volume=true;
            
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
              //cout<<"CC-";
              //CC_num++;
              //cout<<"charged current of neutrinos happened \n";
              // Obtain Bjorken y
              if(part_type==16)
              {
                if(anti==1) tCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                if(anti==-1) tCCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                tau_x=pos[0];
                tau_y=pos[1];
                tau_z=pos[2];
              }
              if(part_type==14||part_type==12)
              {
                if(anti==1)mCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                if(anti==-1)mCCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                if(part_type==12)broken=true;
	      }
              //if(part_type==12)
              //{
              //  broken=true;
                
                //particle becomes e so save in event and break
              //}

              Bjorken_y=finalstatecc[1];//INEASTIVITY?
            
              // Set the tau lepton energy from the sampled Bjorken y.
              double initial_energy=part_energy;
              part_energy=(1.-Bjorken_y)*part_energy;
              double shower_energy=initial_energy-part_energy;

              // Increment the cc interaction counter in the event structure.
              //event.ncc++;
              //add config.save_nu_ev
              if(config.save_nu_events&&in_volume(pos[0],pos[1],pos[2]))
              {
              //if(config.save_sto&&lep.Ei==0)lep.fill_initial(pos[0],pos[1],pos[2],part_energy,anti,part_type);
              //temp_lep.xi=pos[0];
              //temp_lep.yi=pos[1];
              //temp_lep.zi=pos[2];
              //temp_lep.Ei=part_energy;
              //if(part_type==12)Bjorken_y=1;
              got_event=true;
              event_count++;
              //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
              outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<< initial_energy<<","<<Bjorken_y<<","
              <<part_type*anti<<","<<0<<","<<0<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<endl;
              }
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
              
              if(part_type==12)part_type=11;
              if(part_type==14) part_type=13;//nu_mu -> mu
              if(part_type==16) part_type=15;//nu_tau -> tau
              //if(in_volume(pos[0],pos[1],pos[2])&&config.save_sto) lep.p_type=part_type;//temp_lep.p_type=part_type;

              // get the density at the current location before jumping to the tau lepton part of the loop
              //dens = get_dens_from_coords(pos);//change to Lmax2 dfor icecube
              change=true;//signify that particle type changed so simulation of decay delayed until next distance loop iteration
              //ADD EVENT SAVING!!
              CC_num++;
              if(broken==true)break;
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
              
              if(part_type==16)
              {
                if(anti==1) tNCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                if(anti==-1) tNCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
               
              }
              if(part_type==14 || part_type==12)
              {
                if(anti==1) mNCFinalData->ThrowFinal(log10(part_energy),finalstatenc);
                if(anti==-1) mNCBarFinalData->ThrowFinal(log10(part_energy),finalstatenc);
              }
              
              Bjorken_y=finalstatenc[1]; 
              
              // Set the neutrino energy from the sampled Bjorken y.
              double initial_energy=part_energy;
              part_energy=(1.-Bjorken_y)*part_energy;
              double shower_energy=initial_energy-part_energy;


              

              generation++;
              //add in config.save_nu_ev
              if(config.save_nu_events&&in_volume(pos[0],pos[1],pos[2]))
              {
              //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
              outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_energy<<","<<Bjorken_y<<","
                  <<part_type*anti<<","<<1<<","<<0<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<endl;
              event_count++;
              got_event=true;
              }
              NC_num++;

              change=true;

            }
            else if (GRhappens)
            {
              //total_GR++;
              //GR_num++;
              double initial_energy=part_energy;
              //cout<<"GR-";
              if(part_type!=12 || anti!=-1) 
              {
                cout<<CCratio<<" "<<NCratio<<" "<<GRratio<<endl;
                cout<<"CC happens = "<<CChappens<<" NChappens = "<<NChappens<<" GR happens = "<<GRhappens<<endl;
                cout<<"GR without a valid particle"<<endl;
                continue;
              }
              
              double react_rand=((double) rand() / (double)(RAND_MAX));
              //cout<<react_rand<<endl;
              int temp_channel=-1;
              if(react_rand<0.676)
              { 
                temp_channel=0;
                //cout<<"W+ decayed to quarks"<<endl;
                //add in config.save_nu_ev
                if(config.save_nu_events&&in_volume(pos[0],pos[1],pos[2]))
                {
                //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
                outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_energy<<","<<1<<","
                    <<part_type*anti<<","<<2<<","<<temp_channel<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<endl;
                
                event_count++;
                got_event=true;
                }
		GR_num++;
                broken=true;
                break;
              }
              else
              {
                
                //cout<<"W+ decayed to leptons"<<endl;
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
                  temp_channel=1;
                  
                }
                else if(lep_rand<(2./3.) &&lep_rand>=(1./3.))
                {//mu made
                  part_type=14;
                  lepton_mass=mmuon;
                  lepton_type=13;
                  temp_channel=2;
                  
                }
                else if(lep_rand>=(2./3.))
                {//tau made
                  part_type=16;
                  lepton_mass=mtau;
                  lepton_type=15;
                  temp_channel=3;
                  
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
                //add in config.save_nu_ev
                if(config.save_nu_events&&in_volume(pos[0],pos[1],pos[2]))
                {
                outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_E<<","<<gr_inel<<","
                    <<part_type*anti<<","<<2<<","<<temp_channel<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<endl;
                event_count++;
                //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
                got_event=true;
                }
               
                

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
                  if(config.save_sto)particle_data.start_in_volume.push(in_volume(pos[0],pos[1],pos[2]));
                }
                
              }
              change=true;//signify that particle type changed so simulation of decay delayed until next distance loop iteration
              GR_num++;
            }
          }// end of 'if(Ldist<Lmax)'  (particle still inside Earth)

        }// end of 'if (tag == 0,1,2)' (particle is a nu_e,nu_m,nu_t)
        
        //=========================
        // Particle is a tau lepton or muon lepton
        //=========================
        if((part_type==11||part_type==13||part_type==15)&&change==false&&broken==false) //change particle types
        {


          //bool cont_sto=false;//temp var for which prop moldule to use
          double frac_loss=0;
          double sampled_decay_length=0;
          double random_num=0;

          //if(!in_volume(pos[0],pos[1],pos[2]))
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
          random_num=(double)rand()/(double)RAND_MAX;
          sampled_decay_length=decay_length((double)rand()/(double)RAND_MAX,part_energy,part_type); //cm

          //get the interaction info
          if(config.use_sto_inst_of_cont) //true is sto
          {
            sto.set_val(part_energy,part_type,dens);
            dL=sto.get_interaction_length();
            frac_loss=sto.get_sampled_energy();
          }
          else //false is cont
          {
            cont.set_values(part_energy,dens,part_type,0);
            dL=cont.get_interaction_length();
            frac_loss=cont.get_energy_loss()/part_energy;
          }



          //if(dL>1E6){cout<<"dens"<<dens<<"   buggy particle from a super long step with step "<<dL<<"and en loss "<<frac_loss<<endl;}
          // Estimate step length based on Energy, dE/dx, local density, and fraction.
          //dL=(part_energy/(dens*elost(part_energy, dens, ELOSSmode,part_type)))*frac;
          
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
          
          // Calculate the probability that the tau lepton will decay along the step dL.
          //dPdes=1.-exp(-dL*dPdesdx(part_energy,part_type));
          
          //rndm=((double) rand() / (double)(RAND_MAX));
          
          //update particle position
          //if(rndm<dPdes)cout<<"tau to decay here "<<part_pos<<endl;
          //cout<<part_energy<<","<<frac_loss<<","<<part_energy*frac_loss*pow(10,9)<<","<<config.energy_threshold<<endl;
          //cout<<frac_loss<<endl;
          
          if(config.save_sto_events&&in_volume(pos[0],pos[1],pos[2])&&frac_loss*part_energy>Elim)//save the deposition
          {
            //output this
            //cout<<"stuff here"<<endl;
            int temp_show=0;
            if(sto.sto_type==0)temp_show=1;
            if(sto.sto_type==1)temp_show=1;//sto_type==0 brem, sto_type==1 pp, sto_type==2 pn
            if(sto.sto_type==2)temp_show=0; 
            outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<part_energy<<","<<frac_loss<<","
              <<part_type*anti<<","<<sto.sto_type+4<<","<<temp_show<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<sto_num<<endl;
            //outLep<<log10(part_energy)<<","<<frac_loss<<","<<log10(part_energy*frac_loss)<<","<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<sto.sto_type<<","<<i<<","<<num_count<<","<<part_type<<","<<sto_num<<endl;
            sto_num++;
            //log_initial_energy,frac_loss,log_dep_energy,x,y,z,vx,vy,vz,int type, traj id, part type
          }
          //dont flatly update. instead depends on decay length
          //part_pos=part_pos+dL;
          //pos[0]=pos[0]+dL*x_step;
          //pos[1]=pos[1]+dL*y_step;
          //pos[2]=pos[2]+dL*z_step;

          //part_energy=(1.-frac_loss)*part_energy;

          //if(rndm<dPdes)cout<<"tau decays here "<<part_pos<<endl;

          if(in_volume(pos[0],pos[1],pos[2])&&!checked_in)
          {
            taus_passed_through++;
            checked_in=true;
          }
          if(in_volume(pos[0],pos[1],pos[2])&& !in_vol)
          {
            entered_volume=true;
            
          }
        
          // Calculate a random number used for determining whether the tau decays or not in this step.
          

          //part_energy = part_energy-dL*dens*elost(part_energy, dens, ELOSSmode,part_type);
          // Determine whether the lepton decays or not.
          
          //if(config.save_sto&&in_volume(pos[0],pos[1],pos[2])&&lep.Ei==0)
          //{
          //  lep.fill_initial(pos[0],pos[1],pos[2],part_energy,anti,part_type);

          //}

          if(sampled_decay_length>dL)
          {
            //==============================
            // The tau lepton does NOT decay
            //=============================
            //cout<<"no decay"<<endl;
            part_pos=part_pos+dL;
            pos[0]=pos[0]+dL*x_step;
            pos[1]=pos[1]+dL*y_step;
            pos[2]=pos[2]+dL*z_step;

            part_energy=(1.-frac_loss)*part_energy;

            if(part_pos+dL > maxL) 
            {
              dL=maxL-part_pos;//change tolmax1 for icecube
              //cout<<"tau left before decaying "<<endl;
            }
            // Calculate the traversed grammage
            traversed_grammage+=dL*dens;
            
            //cout<<"ND-";
          }
          else
          {
            //cout<<"D-";
            save_lep=true;
            if(in_volume(pos[0],pos[1],pos[2]))
            {
              decay_num++;
            }

            part_pos=part_pos+sampled_decay_length;
            pos[0]=pos[0]+sampled_decay_length*x_step;
            pos[1]=pos[1]+sampled_decay_length*y_step;
            pos[2]=pos[2]+sampled_decay_length*z_step;

            if(part_pos+sampled_decay_length > maxL) 
            {
              sampled_decay_length=maxL-part_pos;//change tolmax1 for icecube
              //cout<<"tau left before decaying "<<endl;
            }
            // Calculate the traversed grammage
            traversed_grammage+=sampled_decay_length*dens;
          
            //=======================
            // The lepton decays
            //=======================
            //if(config.save_sto&&in_volume(pos[0],pos[1],pos[2]))
            //{
            //lep.fill_end(pos[0],pos[1],pos[2],part_energy);
            //lep.save_lep(&outLep,0,i,num_count);
            
            //cout<<"charged in volume - should be a saved lep"<<endl;
            //}
            //if(config.save_sto)lep.empty_class();
            //dc_num++;
            //cout<<"tau path length before decaying is "<<tau_path<<endl;
            //cout<<"tau energy at decay "<<part_energy<<endl;
            //cout<<"dL is "<<dL<<endl;
            //cout<<"dpdes is "<<dPdes<<endl;
            //cout<<"rndm is "<<rndm<<endl;
            
            //if(part_type==6) num_tau_decays++;
            //if(part_type==5) num_muon_decays++;
            
              //temp_lep.xf=pos[0];
              //temp_lep.yf=pos[1];
              //temp_lep.zf=pos[2];
              //temp_lep.Ef=part_energy;
              //save_lep=true;
            
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
            int is_had_or_em=-1;
            if(part_type==15)
            { 
              //cout<<"tau decay"<<endl;

              if(config.conversion)
              {
                for(int j=1;j<6;j++)
                {
                  if(reaction_data.tau_type[reaction_index][j]!=0) //
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
                  if(is_had_or_em==-1 && reaction_types[j]==0) is_had_or_em=0;//for EM
                  if(is_had_or_em==-1 && reaction_types[j]==11) is_had_or_em=1;//for EM
                  if(is_had_or_em==-1 && reaction_types[j]==13) is_had_or_em=2;//for EM
                }
              //ADD IF HAD OR EM CHANNEL, EITHER 12,13,14,15 or 0
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
                  if(is_had_or_em==-1 && reaction_types[j]==0) is_had_or_em=0;//for HAD
                  if(is_had_or_em==-1 && reaction_types[j]==11) is_had_or_em=1;//for EM
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
            //add in config.save_dec
            if(config.save_dec_events&&in_volume(pos[0],pos[1],pos[2]))
            {
            outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_energy<<","<<frac_energy_dumped<<","
                <<initial_particle*anti<<","<<3<<","<<is_had_or_em<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<endl;
            event_count++;
            //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
            got_event=true;
            }
            dc_num++;
            
            if(config.conversion)
            {
              for( int j=1;j<6;j++)
              {
                  
                if((reaction_energies[j]>Elim)&&(part_pos<Lmax)&&(reaction_types[j]!=0)) //particles inside earth and above threshold are stacked, change to lam2 for icecube
                {
                  //cout<<"pushing new particle to stack \n";
                  //produced_muons+=1;
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
                  if(config.save_sto)particle_data.start_in_volume.push(in_volume(pos[0],pos[1],pos[2]));
                }
              }
            }
            if(!config.regen) { broken=true; break;}
          }
        }// end of 'if (tag == 1,2,3)' (particle is lepton)
          


       

        //time to pop all the extra particles back in the stack to be looped over after
        
        
        if(in_vol&&!in_volume(pos[0],pos[1],pos[2]))
          {
            //cout<<"left volume"<<endl;
            left_volume=true;
            
          } 
       
        if(part_energy<Elim){
          //if(part_type==3) tau_neutrinos_below++;
          //else if (part_type==2) muon_neutrinos_below++;
          //else if (part_type==6) taus_below++;
          //else if (part_type==5) muons_below++;
          
          broken=true;
          break;
          }
      } // ends loop 'for(Ldist=0.;Ldist<Lmax;)'
      
      //=================================================
      // Write energy of emerging tau to output text file
      //=================================================
      //cout<<"save- type: "<<part_type<<" .energy: "<<part_energy<<".pos: "<<part_pos<<". gen: "<< generation<<endl;
    

      for(int j=0; j<5;j++)
      {
        //cout<<type_to_save[i]<<","<<part_type<<endl;
        //cout<<"part energy is "<<part_energy<<"with threshold "<<Elim<<endl;
        if(config.save_emerging&&(part_type==type_to_save[j])&&(part_energy>Elim)&&broken==false&&has_been_saved==false)//changge to just poarticle tyoe if saving all particls even below thrwhpold
          {//add &&false after testing influx and outflux files
            out_leptons++;
            //cout<<"got one!"<<endl;
	          outEnergies <<part_type<<" "<<anti<<" "<<NC_num << " " << CC_num << " " << GR_num<<" "<<
            dc_num << " " <<generation <<" "<<num_count<< " " << log10(part_energy)+9 << " " << log10(Energy_GeV)+9<<" "<<part_pos<<" "<<tau_x<<" "<<tau_y<<" "<<tau_z<<"\n";
            has_been_saved=true;
            
            if(out_leptons==500)outEnergies<<num_count<<" initial neutrinos of type "<<config.starting_type<<" at energy "<<Energy_GeV*pow(10,9);
          }
        
      }//if the main particle isn't saved then it is forgetten falling below thrshold
      
    } while(!particle_data.part_type.empty()); // end of loop of stack
    }
  delete cumulative_grammage;
  delete grammage_distance;
  
  } // end of loop for initial particles
  

  printf("taus decayed: %i, taus passing through (could have decayed): %i",decay_num,taus_passed_through);
  outEnergies << "END" << endl; // write END in the last line of the text output file. 
  outEnergies.flush();
  outEvents <<"END"<<endl;
  outEvents.flush();
  //out_gram<<"END"<<endl;
  //out_gram.flush();

  //outLep<<"END"<<endl;
  //outLep.flush();
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
  
  return 0;
  
}    

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// ===================================================
// Several functions
// ===================================================

double decay_length(double P, double E, int particle_type)
{
    double mass=0;
    double lifetime=0;
    if(particle_type==15)
    {
        mass=mtau;
        lifetime=tau_lifetime;
    }
    else if(particle_type==13)
    {
        mass=mmuon;
        lifetime=muon_lifetime;
    }
    else
    {
        cout<<"wrong particle"<<endl;
        exit(-1);
    }
    if(1-P<0)
    {
      cout<<"log bug"<<endl;
      exit(-1);
    } 

    double decay_length=-log(1-P)*E*c_light*lifetime/mass*1e2;
    return decay_length; //cm
}


void set_points_from_angle(input_file *input,double angle)
{
  double x_step=cos((angle-90)*PI/180);
  double y_step=sin((angle-90)*PI/180);
  
  double start_point[2];
  double end_point[2];
  double depth=0;
  double small_r=0;
  end_point[0]=0;
  end_point[1]=R0-depth+small_r*y_step;
  cout<<"end point "<<end_point[0]<<","<<end_point[1]<<endl;
  double point_slope=0;
  double y_pos=0;
  double y_neg=0;

  if(angle!=180.0||angle!=0)
  {
    point_slope=-tan((angle-90)*PI/180);
    double det=4*(R0-depth)*(R0-depth)-4*(1+point_slope*point_slope)*((R0-depth)*(R0-depth)-R02*point_slope*point_slope);
    double denom=2*(1+point_slope*point_slope);
    if(det<0)
    {
      cout<<"invalid, breaking code"<<endl;
      //return -1;
    }
    y_pos=(2*(R0-depth)+sqrt(det))/denom;
    y_neg=(2*(R0-depth)-sqrt(det))/denom;
  }
  if(angle>=90.0)start_point[1]=y_neg;
  if(angle<90.0)start_point[1]=y_pos;
  
  start_point[0]=sqrt(R02-start_point[1]*start_point[1]);
  cout<<"start point "<<start_point[0]<<","<<start_point[1]<<endl;
  if(angle==180.0)
  {
    start_point[1]=-R0;
    start_point[0]=0;
  }
  if(angle==0.)
  {
    start_point[1]=R0;
    start_point[0]=0;
  }
  input[0].eex=0;
  input[0].eey=0;
  input[0].eez=R0;

  input[0].xe=start_point[0];
  input[0].ye=0;
  input[0].ze=start_point[1];

  input[0].xi=start_point[0];
  input[0].yi=0;
  input[0].zi=start_point[1];

  input[0].xf=0;
  input[0].yf=0;
  input[0].zf=R0;
  }
bool in_volume(double x,double y,double z) //(det* detec, double x,double y,double z)
{
  double rad=sqrt(x*x+y*y);
  if(rad<=15*pow(10,5)&&(z<R0)&&(z>R0-2.8*pow(10,5)))
  {
    return true;
  }
  else return false;
  /* for when transitioning to detec config file
  if(detec.inner_rad==0)
  {
   double sp_R=sqrt(x*x+y*y+z*z);
   if(sp_R<detec.inner_sphere) return true;
   else return false;
  }
  if(detec.inner_phere==0)
  {
     double rad=sqrt(x*x+y*y);
    if(rad<=15*pow(10,5)&&(z<R0)&&(z>R0-2.8*pow(10,5)))
    {
      return true;
    }
    else return false;
  }
  */
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
    if(row_ind+1>file_length) 
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
  
  return nameEnergies;
}
string make_event_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp, int p_type)
{
  string type_temp="";
  switch (p_type){
    case 12: type_temp="nue"; break;
    case 14: type_temp="numu"; break;
    case 16: type_temp="nutau"; break;
    case -12: type_temp="anti_nue"; break;
    case -14: type_temp="anti_numu"; break;
    case -16: type_temp="anti_nutau"; break;
  }
  string nameEvents="";
  nameEvents+=config.data_dir;
  nameEvents+="/events/";
  nameEvents+=type_temp;

  nameEvents+="_events_";
  nameEvents+=es_temp;
  //nameEvents+="_";
  //nameEvents+=angs_temp;
  nameEvents+=".dat";

  return nameEvents;

}

string make_lepton_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp, int p_type)
{
  string type_temp="";
  switch (p_type){
    case 12: type_temp="nue"; break;
    case 14: type_temp="numu"; break;
    case 16: type_temp="nutau"; break;
    case -12: type_temp="anti_nue"; break;
    case -14: type_temp="anti_numu"; break;
    case -16: type_temp="anti_nutau"; break;
  }
  string nameLeptons="";
  nameLeptons+=config.data_dir;
  nameLeptons+="/leptons/";
  nameLeptons+=type_temp;

  nameLeptons+="_leptons_";
  nameLeptons+=es_temp;
  //nameLeptons+="_";
  //nameLeptons+=angs_temp;
  nameLeptons+=".dat";

  return nameLeptons;

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
     else if ((int)line.find("n_throws")!=-1) sin >>config.n_throws;
     else if ((int)line.find("n_traj")!=-1) sin >>config.n_traj;
     else if ((int)line.find("save_sto")!=-1) sin >>config.save_sto;
     else if ((int)line.find("save_nu_events")!=-1)sin>>config.save_nu_events;
     else if ((int)line.find("save_sto_events")!=-1)sin>>config.save_sto_events;
     else if ((int)line.find("save_dec_events")!=-1)sin>>config.save_dec_events;
     else if ((int)line.find("save_emerging")!=-1)sin>>config.save_emerging;
     else if ((int)line.find("use_sto_inst_of_cont")!=-1)sin>>config.use_sto_inst_of_cont;

  }
}

void load_geo() //could be useful later on
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
     else if ((int)line.find("outer_rad")!=-1) sin>>det.outer_rad;
     else if ((int)line.find("outer_depth")!=-1) sin>>det.outer_depth;
  }
  if(det.outer_rad==0) det.traj_weights=(2*det.outer_rad*PI*PI + det.outer_depth*det.outer_rad*2*PI)/R02;
  if(det.outer_sphere==0) det.traj_weights=2*PI*det.outer_sphere*det.outer_sphere/R02;
  
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
/*
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
*/
/*
void generate_events(det * detec,int num_traj)
{
  string output_file="in_files/sim_generated_traj.csv";
  ofstream out_traj(output_file.c_str());
  double x,y,z,u,v,w;

  if(detec.outer_rad==0)
  {

  }
  if(detec.outer_sphere==0)
  {
    for (int di=0;di<num_traj;di++)
    {
        z=(double)rand()/(double)RAND_MAX*detec.inner_depth;
        x=(double)rand()/(double)RAND_MAX*detec.inner_rad;
        y=(double)rand()/(double)RAND_MAX*sqrt(detec.inner_rad*detec.inner_rad-x*x);
        u=(double)rand()/(double)RAND_MAX*1;
        v=(double)rand()/(double)RAND_MAX*sqrt(1-u*u);
        w=(double)rand()/(double)RAND_MAX*sqrt(1-u*u-v*v);

        //step forward to find earth entrance and stepo backward to find earthexit
    }

  
  }


  out_traj.flush();
  out_traj.close();


}
*/
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
  mkdir((dirs+"/leptons").c_str(),0755);
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



