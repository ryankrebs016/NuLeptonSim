//############################################################################# 
// Propagation of neutrinos in Earth - flux of emerging leptons
//
// Adapted from code used to obtain probability of emerging tau leptons in 
// calculations of exposure to Earth-skimming neutrinos in Auger
//############################################################################# 
// Processes:
// - CC or NC interaction of nu_tau and nu_mu in Earth (variable density along chord) 
// - Production of tau and muon leptons with sampling of (1-y) where E_tau=(1-y)*E_nu_tau 
// - Propagation of leptons (including energy loss)
// - Tau decay sampled from pythia look up tables for resulting neutrino and lepton-neutrino pairs
// - Reinteraction of nu_tau produced in tau decay 
// - Reinteraction of nu_tau produced in nu_tau NC interaction
//----------------------------------------------------------------------------- 
// Several models of neutrino cross-section & lepon energy loss can be chosen
//----------------------------------------------------------------------------- 
// Energies (GeV), unless otherwise specified.
//----------------------------------------------------------------------------- 
// Assume all particles are not anti-matter
// No electron generation particles are simulated yet
// Config file can be changed to consider regeneration and decays
// Config files stores settings used in the code such as the starting neutrino type,
// bool to use energy distribution, what angles to simulate over
//-----------------------------------------------------------------------------
// calling the code is the form ./Simu_elost 1E+20 95.0 1E+4 0 0 4.0 0.92 test ./
// if energy dist is used or angles are simulated over then the two values for them is ignored

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

#include "Table.hh"
#include "Earth.hh"
#include "Constantes.hh"

using namespace std;



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

typedef struct {

  stack<int> part_type;       // 0-5
  stack<double> part_energy;  // E_min<part_energy<E_init
  stack<double> part_pos;     // position inside earth
  stack<bool> anti;           // true if anti particle
  stack<int> generation;      // number of interactions before this particle was produced
  stack<int> NC_num;          // number of NC preceding creation of this particle
  stack<int> CC_num;          // number of CC preceding creation of this particle
  stack<int> dc_num;          // number of decays it took to get to the current particle

} particle_info_def;

typedef struct {

  int tau_type[10000][6];       // hold particle types created in tau decay
  double tau_energy[10000][6];  // hold energy of created particle  in tau decay
  int mu_type[10000][6];        // hold particle types created in muon decay
  double mu_energy[10000][6];   // hold energy of created particles in muon decay
    
} reaction_tables_def;

typedef struct {
  int starting_type;          // starting type of neutrinos
  bool sim_angles;            // true or false to simulate all angles in the range
  double low_angle;           // smallest angle in the simulation range
  double high_angle;          // largest angle in the simulation range
  bool regen;                 // true or false for lepton regeneration through decays
  bool conversion;            // true or false to include produced particles from W +/- decays
  bool energy_distribution;   // trueo or false to use energy distribution

}config_init;   // data struct to hold the value read from config file

MYEVT_DEF event;
particle_info_def particle_data;
reaction_tables_def reaction_data;
config_init config;



// initialize reaction from pythia table and convert pythia type to tag code used in code
void initialize_reaction(int tau_type[][6], int mu_type[][6],double tau_ene[][6],double mu_ene[][6]);

// converts pythia tags to tags in this code
int convert_types(int pythia_type);

//load values from cpnfig file
void load_config();

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
// Tau and muon neutrino cross sections: CC and NC
double dsigCC(double l1, int CCmode, int type, int AntiNu);	
double dsigNC(double l1, int CCmode, int type, int AntiNu);

// -------------------------------------------------
// Local density as a function of zenith angle
double earthdens( double *x, double *par); 

// -------------------------------------------------
// Average density as a function of zenith angle
double mean( double *x, double *par);
double mean_dens_chord(double theta);

//---------------------------------------------------------------
// Declare event structure - contains information about event: 
// nu_tau propagating and creating new nu_tau in NC or tau in CC
// where tau can decay and produce a new nu_tau.
//---------------------------------------------------------------
//---------------------------------------------------------------

// Initialize Earth class
// The arguments are water thickness and density. 
// They are initialized to bare rock here but it is re-initialized below.
Earth *terra = new Earth(0.0, 2.6); 


//#############################################################
// Main code
//#############################################################
int main(int argc, char **argv)
{
  load_config(); // call function to load config
 
  // loads reaction data from pythia table
  for(int i=0;i<10000;i++){for(int j=0;j<6;j++){reaction_data.tau_type[i][j]=reaction_data.mu_type[i][j]=-1;reaction_data.tau_energy[i][j]=reaction_data.mu_energy[i][j]=0.0;};}  
  initialize_reaction(reaction_data.tau_type,reaction_data.mu_type,reaction_data.tau_energy,reaction_data.mu_energy);
  
  //bool sim_all_angles=false;
  //bool consider_conversion= false; // generation conversion - electrons are still neglected bc of energy loss snd cross sectional models are needed.
  int n1=(95-config.low_angle)/0.1; // changing second value of 90 to 95 so it can run w/ icecube
  int n2=(config.high_angle-95)/1;
  //int size=n1+n2;
  //double ang[n1+n2];
  stack <float> angles;
  //int starting_type=2; // 0 = Ve , 1= Vm , 2= vT, 3=e, 4=m, 5=T
  int type_to_save[4]={config.starting_type,-1,-1,-1};
 
  if(config.regen)
  {
    type_to_save[1]=config.starting_type+3;
  }

  if(config.conversion)
  {
    if(config.starting_type==1) {type_to_save[2]=2; type_to_save[3]=5;}
    if(config.starting_type==2){type_to_save[2]=1; type_to_save[3]=4;}
  }

  if(config.sim_angles)
  {
    for (int i=n1;i>0;i--) angles.push(95.-i*0.1); 
    for (int i=0;i<=n2;i++) angles.push(95.+i);  
  }
  else
  {
    angles.push(atof(argv[2]));
  }
  for(int i=0;i<4;i++) cout<<type_to_save[i]<<" ";
  cout<<endl;

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
  
  cout << "Tau Data Table Initialized: " << InitTau << endl << endl;
  cout << "Moun Data Table Initialized: " << InitMuon << endl << endl;

  bool useEnergyDistribution = false;
  if (atof(argv[1]) == 0) {
    useEnergyDistribution = true;
    cout << "Will throw uniformly random x neutrinos energy between log10(E_nu/eV) = 15 and  log21(E_nu/eV)" << endl;
  }
  // The finalstate array gets filled by sampling of TauData
  double finalstate[6];      // 0=nu_tau, 1=nu_mu, 2=nu_e, 3=hadron, 4=muon, 5=electron
  
  // Initialize Random number generator.
  struct timeval time_struct;
  gettimeofday(&time_struct,NULL);
  srand((time_struct.tv_sec * 1000) + (time_struct.tv_usec / 1000));
  //cout << "Check for random number uniqueness." << endl;
  //for (int ii = 0; ii<3; ii++)
  //{
  //    cout << "Random Test " << ((double) rand() / (double)(RAND_MAX)) << endl;
  //}
  
  while(!angles.empty())  //loops over all angles in the stack, should be a single angle or a stack of angles over a range
  {
    double angle=angles.top();
    angles.pop();
  
    // Cut in energy below which particles are no longer propagated
    double Elim_eV=1.e14;       // eV
    double Elim=Elim_eV*1.e-9;  // GeV
    
    //-------------------------------------------------
    // Declare several variables
    double rndm;           // random number throughout the code
    
    double Ldist;         // Distance along chord length in Earth
    double Depth1=2.45e5;          // Depth of IceCube bottom layer
    double Depth2=1.45e5;          // Depth of IceCube top layer
    double Lmax1;          // Chord length in Earth to bottom layer
    double Lmax2;          // Chord length in Earth to top layer
    double Lmax;
    double dL=0.;	 	  // Propagation step along chord length. Initialized to zero but it varies for each step.
    double frac=1e-3;	  // Fraction of energy lost by tau in each step. This value stays constant
    double dPdes;         // Probability of tau decay in dL
    double traversed_grammage; // Track the grammage traversed upon entering the Earth.
    double Energy_GeV;	  // Particle energy
    double Bjorken_y;     // Bjorken y value for interactions CC & NC
    
    int brkcnt=0, brkcnt2=0; // These are used to flag particles that have fallen below the energy threshold (Elim)
    int prop_mode=0;	  // Used in tau propagation
    
    int AntiNu=0;

    double finalstatecc[2],finalstatenc[2]; // will contain Bjorken (1-y) and y
    
    float dens; // local density along trajectory
    
  
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    // Change following parameters depending on application
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    // Number of neutrinos simulated at each zenith angle
    int tot_evt=1e7;
    if(argc>3){
      tot_evt=(int)atof(argv[3]);
    }
    
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    cout << "======================================" << endl;
    cout << "Number of neutrinos simulated = " << tot_evt << endl;
    cout << "Energy = " << atof(argv[1]) << " eV" << endl;
    cout << "Threshold energy = " << Elim*1.e9 << " eV" << endl;
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    //for rounding energies and angle
    string e_temp=to_string(log10((double)atof(argv[1])));
    string es_temp="";
    string ang_temp=to_string(angle);
    string angs_temp="";
    int count=0;
    for(int i =0;i<4;i++)es_temp+=e_temp[i];
    
    if(angle>=100) count=5;
    else count =4;
    for(int i=0;i<count;i++) angs_temp+=ang_temp[i];
      

    
  

    //-------------------------------------------------
    // Output file names using input arguments
    string nameEnergies="";
    /*
    if(argc>=9){
      nameEnergies+=argv[9];
    }
    */
    nameEnergies+="data/";
    nameEnergies+="particles_";
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
    ofstream outEnergies(nameEnergies.c_str());
    outEnergies << "10^"<< log10(tot_evt)<<" initial neutrinos of type "<<config.starting_type<<" at energy "<<argv[1]<<". \nKey: type,NC,CC,DC,Gen,Energy.\n";
    
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
      cout << "Outer Layer Thickness " << terra->depth_new_layer   << " km" << endl;
      cout << "Outer Layer Density   " << terra->dens_new_layer    << " g/cm^3" << endl;
    }
    
    // Get zenith angle from input argument
    double refTheta=angle;
    cout << "Theta " << refTheta << " deg" << endl;
    
    double min_depth = R0-R0*sin(PI*(1.-refTheta/180.));
    if(Depth2>min_depth) 
    {
      cout << "Angle too steep to be measured by IceCube" << endl;
      break;
    }

    Lmax1=R0*cos(PI*(1.-refTheta/180.))+sqrt(pow(R0*cos(PI*(1.-refTheta/180.)), 2.0)-(2*R0*Depth1-pow(Depth1,2.0))); // chord length inside Earth in cm
    Lmax2=R0*cos(PI*(1.-refTheta/180.))+sqrt(pow(R0*cos(PI*(1.-refTheta/180.)), 2.0)-(2*R0*Depth2-pow(Depth2,2.0))); // chord length inside Earth in cm
    // Earth's chord (km) at thetaRef (deg) (R0 is radius of Earth in cm)
    cout << "Chord length Bottom " << Lmax1/pow(10,5) << "km" << endl;
    cout << "Chord length Top " << Lmax2/pow(10,5) << " km" << endl;
    
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
    
    if(Lmax2<=0.) Lmax2=0.;                 // set to zero  if the value is negative.
    
    printf ("Lmax %1.3e cm, %1.2f deg \n", Lmax2, refTheta);
    
    // create look-up datable of distance as a function of grammage
    double sum_grammage = 0.;
    for (int ii=1; ii<=1000000; ii++)
    {
    double dx = Lmax2/1000000;
    double x_val = Lmax2 * double(ii) / 1000000.;
    sum_grammage +=  dx*earthdens(&x_val, &Lmax2);
    }
    printf("sum_grammage %1.5e g/cm^2\n",sum_grammage);
    printf("mean_dens %1.2e\n\n",sum_grammage/Lmax2);
    
    double* cumulative_grammage = new double[1000000];
    double* grammage_distance = new double[1000000]; // in cm
    
    double d_grammage = sum_grammage/1000000.; // g/cm^2
    
    cumulative_grammage[0] = 0.;
    grammage_distance[0]   = 0.;
    for (int ii=1; ii<=1000000; ii++)
    {
    double x_val = grammage_distance[ii-1];
    double dx    = d_grammage/earthdens(&x_val, &Lmax2);
    cumulative_grammage[ii] = cumulative_grammage[ii-1] + dx*earthdens(&x_val, &Lmax2);
    grammage_distance[ii] = x_val + dx;
    //printf("*** ii %d %1.5f %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
    //if(ii%100000 ==0) printf("ii %d %1.2e %1.5f\n",ii, grammage_distance[ii], cumulative_grammage[ii]);
    }
    // start with one particle then it will loop over stack and iterate to the next original particle
    for( int i=0;i<tot_evt;i++)
    {
      //cout<<endl<<"original particle: "<<i<<endl;
      Energy_GeV = atof(argv[1])*pow(10,-9); // Get nu_tau energy from input argument
      if (config.energy_distribution) Energy_GeV =  pow(10,6 + (6 * (double)rand()/RAND_MAX));
      //cout<<"initial energy GeV is "<<Energy_GeV<<endl;
      //cout<<"threshold energy in GeV is "<<Elim<<endl;
      int part_type=config.starting_type;
      double part_energy=Energy_GeV;
      double part_pos=0;
      bool anti=0;
      int generation=0;
      int NC_num=0;
      int CC_num=0;
      int dc_num=0;
      /*
      particle_data.part_type.push(part_type);
      particle_data.part_energy.push(part_energy);
      particle_data.part_pos.push(part_pos);
      particle_data.anti.push(anti);
      particle_data.generation.push(generation);
      particle_data.NC_num.push(NC_num);
      particle_data.CC_num.push(CC_num);
      particle_data.dc_num.push(dc_num);
      */
      //cout<<"Initial- type: "<<part_type<<" .energy: "<<part_energy<<"GeV. pos: "<<part_pos<<"cm . gen: "<< generation<<endl;
    
      int loop_num=0;
      //======================= Loop over stacks until empty
      do 
      { 
        if(loop_num!=0)cout<<"stack loop number:"<<loop_num<<endl;
        loop_num++;

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
          //cout<<part_type<<endl;
          //======================================================
          particle_data.part_type.pop();
          particle_data.part_energy.pop();
          particle_data.part_pos.pop();
          particle_data.anti.pop();
          particle_data.generation.pop();
          particle_data.NC_num.pop(); 
          particle_data.CC_num.pop();
          particle_data.dc_num.pop();
        }
        int reaction_types[6]={-1,-1,-1,-1,-1,-1};
        double reaction_energies[6]={-1,-1,-1,-1,-1,-1};
        int loops=0;
        if(part_type==0||part_type==3) continue; //ignores particles below threshold or of electron generation
        if(part_energy<Elim) continue;
        //create holding arrays for reactions
        

        // Flag to see how fast code is running
        //if(!((float)i/100000-(int)(i/100000))) cout<< i << endl;
        
        // Get nu_tau energy from input argument
      
        //    cout << Energy_GeV << endl;
        //if(starting_type=="muon"){
        //  AntiNu=rand()%2;
        //};
        
        /* ------------------event struct is commented out as it is no longer needed
        // Initialize event info
        //event.nevt = i;            // event number
        event.theta = refTheta;      // zenith angle
        event.Lmax1 = Lmax1;         // max. distance in Earth
        event.Lmax2 = Lmax2;         // max. distance in Earth
        event.Estart = Energy_GeV; // initial nu_tau energy
        event.trig = 0;            // no trigger (tau emerging) so far
        event.ncc = 0;             // zero CC interactions so far
        event.nnc = 0;             // zero NC interactions so far
        event.ndk = 0;             // zero tau decays so far
        event.npart = 1;           // only 1 particle so far
        int tag=0;	               // dealing with a nu_tau
        int npart = 0;
        event.Shheight = 0;
        event.Shlong = 0;
        event.E1[npart] = Energy_GeV;
        event.v1[npart] = 0;
        event.E2[npart] = -1;
        event.v2[npart] = -1;
        event.id[npart]=0;
        */
        brkcnt=0;
        prop_mode =0;
        traversed_grammage = 0.; //initialize traversed grammage for each event
        
        //======================= Start propagation along chord of length Lmax
        //printf("Ldist, Lmax %1.2e %1.2e\n", Ldist, Lmax);
        //Ldist=0.;
        
        for(part_pos;!(((Lmax1<=part_pos) && (part_pos<=Lmax2))||part_pos>=Lmax2);)
        {
          loops++;
          //cout <<"distance is "<< part_pos << "and energy is "<<part_energy<<endl;
          // Get the local density for this part of the chord.
          dens = earthdens(&part_pos,&Lmax2);
          //cout<<"In dist loop - type: "<<part_type<<". energy: "<<part_energy<<". pos: "<<part_pos<<". gen: "<< generation<<endl;
          bool change=false; //is needed so in case a NuTau becomes a Tau it won't be seen by the "is tau" section

          //===========================
          // Particle is a neutrino
          //===========================
          if(part_type==0||part_type==1||part_type==2)
            {
        
            // Number of interaction lengths propagated in this step is given by an exponentially distributed random number.
            double num_int_lens=-log((double) rand() / (double)(RAND_MAX)); // Randomly sampled number of interaction lengths.
            
            // Convert the number of interaction lengths to grammage.
            double X_int;
            X_int = num_int_lens /(Navo*(dsigCC(part_energy, CCmode,part_type,anti)+dsigNC(part_energy, CCmode,part_type, anti)));
        
            // The following lines use the grammage_distance lookup table to estimate what position along the trajectory this Xint corresponds to.
            
            // Add X_int to the total grammage traversed by the particle
            traversed_grammage += X_int;
            
            // Initialize the interaction length distance for this step to zero.
            double Lint = 0.;
          
            // If too large, make sure it exits the volume.
            // NOTE: use floats for this condition. Using ints is bad if float > 2^32, then you get negative int.
            if( traversed_grammage/d_grammage + 1. > 1000000.){
              part_pos = Lmax2; // NOTE: 1000000. is the size of the look-up table.
              traversed_grammage = sum_grammage;
            }
            // If contained within the trajectory, linearly interpolate its interaction distance.
        
            if ( floor(traversed_grammage/d_grammage) + 1. < 1000000.) // NOTE: 1000000. is the size of the look-up table.
            {
              // Get the entry in the look-up table corresponding to the traversed grammage
              int ii_grammage = int(traversed_grammage/d_grammage) + 1;
        
              // Linearly interpolate to estimate the distance propagated
              double slope = (grammage_distance[ii_grammage] - grammage_distance[ii_grammage-1])/d_grammage;
            
              double intercept = grammage_distance[ii_grammage] - slope*cumulative_grammage[ii_grammage];
              Lint = slope*traversed_grammage + intercept - part_pos; // keep track of this step's interaction length.
              part_pos = slope*traversed_grammage + intercept ;
            }
          
            
            // Save the interaction distance propagated in event structure
            //event.L0[npart]=Lint;
            
            // if the neutrino interaction is still inside Earth, simulate an NC or CC interaction and check that particle is still above the tracking energy threshold (Elim)
            if(!((Lmax1<=part_pos) && (part_pos<=Lmax2)))
            {
              // Check if it is CC or NC interaction
              bool CCTauhappens=((double) rand() / (double)(RAND_MAX))>=dsigNC(part_energy, CCmode,part_type,anti)/(dsigNC(part_energy, CCmode,part_type,anti)+dsigCC(part_energy, CCmode,part_type,anti));
              bool CCMuhappens=((double) rand() / (double)(RAND_MAX))>=dsigNC(part_energy, CCmode,part_type,anti)/(dsigNC(part_energy, CCmode,part_type,anti)+dsigCC(part_energy, CCmode,part_type,anti));
              if((CCTauhappens&&part_type==2)||(CCMuhappens&&part_type==1))
              {
                //=======================
                // CC interaction occurs (the tracked particle changes from tau neutrino to tau lepton.)
                //=======================
                
                // Save the energy and position in the event structure
                //event.E2[npart]=Energy_GeV;
                //event.v2[npart]=Ldist;

                
                //cout<<"charged current of neutrinos happened \n";
                // Obtain Bjorken y
                if(part_type==2)tCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
      
                if(part_type==1)
                {
                  if(anti==0)mCCFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                  if(anti==1)mCCBarFinalData->ThrowFinal(log10(part_energy),finalstatecc);
                }

                Bjorken_y=finalstatecc[1];
              
                // Set the tau lepton energy from the sampled Bjorken y.
                part_energy=(1.-Bjorken_y)*part_energy;
                
                // Increment the cc interaction counter in the event structure.
                //event.ncc++;
                NC_num++;
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
                if(part_type==1) part_type=4;
                if(part_type==2) part_type=5;
                // get the density at the current location before jumping to the tau lepton part of the loop
                dens = earthdens(&part_pos,&Lmax2);
                change=true;
                if(!config.regen) break;
              }
              else
              {
                //=======================
                // NC interaction occurs (the tracked particle remains a tau neutrino with reduced energy.)
                //=======================
                
                // Save the energy and position in the event structure
                //event.E2[npart]=Energy_GeV;
                //event.v2[npart]=Ldist;
                //cout<<"neutral current of neutrinos happened \n";
                // Obtain Bjorken y
                if(part_type==2)tNCFinalData->ThrowFinal(log10(part_energy),finalstatenc);  
        
                if(part_type==1)
                {
                  if(anti==0) mNCFinalData->ThrowFinal(log10(part_energy),finalstatenc);
                  if(anti==1) mNCBarFinalData->ThrowFinal(log10(part_energy),finalstatenc);
                }
                
                Bjorken_y=finalstatenc[1];
                
                // Set the neutrino energy from the sampled Bjorken y.
                part_energy=(1.-Bjorken_y)*part_energy;
                
                // Increment the cc interaction counter in the event structure.
                //event.nnc++;           // count NC interaction
                NC_num++;
                                      // Increment the particle counter in the event structure
                //event.npart++;
                //npart++;
                generation++;
              
                // Save the new particle energy and distance in the event structure/
                //event.E1[npart]=Energy_GeV;
                //event.v1[npart]=Ldist;
                
                // Keep the particle tag as neutrino. so no change in part_tupe should happen
                //if(tag==0){event.id[npart]=0;}
                //if(tag==2){event.id[npart]=2;}
                
              }
            }// end of 'if(Ldist<Lmax)'  (particle still inside Earth)
            
            // Energy of new tau or nu_tau produced below threshold => stop
            //if(part_energy < Elim){
              
              // Flag particle below threshold
            //  brkcnt=1;
              
              // Count the number of times the particle fell below threshold
            //  brkcnt2++;
              
            //  break;
            //}
          }// end of 'if (tag == 0,1,2)' (particle is a nu_e,nu_m,nu_t)
          
          //=========================
          // Particle is a tau lepton or muon lepton
          //=========================
          if((part_type==3||part_type==4||part_type==5)&&change==false)
          {
            // Estimate step length based on Energy, dE/dx, local density, and fraction.
            
            dL=(part_energy/(dens*elost(part_energy, dens, ELOSSmode,part_type)))*frac;
            //cout << "Ldist " << 1.e-5*Ldist << "  R " << 1.e-5*sqrt(R02 - (Ldist*Lmax)+Ldist*Ldist) << " dens " << dens << endl;
            //cout<<"simulating tau"<<endl;
            // Check if tau leaves the Earth after dL. If it does then adjust last step
            if(part_pos+dL > Lmax1) dL=Lmax1-part_pos;
            
            // Calculate the traversed grammage
            traversed_grammage+=dL*dens;
            
            // Calculate the probability that the tau lepton will decay along the step dL.
            dPdes=1.-exp(-dL*dPdesdx(part_energy,part_type));
            part_pos=part_pos+dL;
            // Calculate a random number used for determining whether the tau decays or not in this step.
            rndm=((double) rand() / (double)(RAND_MAX));
            part_energy = part_energy-dL*dens*elost(part_energy, dens, ELOSSmode,part_type);
            // Determine whether the tau decays or not.
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
              // Advance the propagation distance by the step dL
              
              cout<<"decay"<<endl;
              // Account for the tau lepton energy lost in the step dL
              
            
              // Save the updated particle energy and location of interaction
              //event.E2[npart]=Energy_GeV;
              //event.v2[npart]=Ldist;
            
              // Get the energy of the neutrino produced in the decay
              generation++;
              //if(tag==1) TauData.ThrowFinal(finalstate);
              //if(tag==3) MuonData.ThrowFinal(finalstate);
              //Energy_GeV=finalstate[0]*Energy_GeV;

              int reaction_index=(double)rand()/(double)RAND_MAX*10000;
              //cout<<"reaction index is:"<<reaction_index<<"\n\n\n";
              //if tau
              if(part_type==5)
              { 
  
                if(config.conversion)
                {
                  for(int i=1;i<6;i++)
                  {
                    if(reaction_data.tau_type[reaction_index][i]!=-1) 
                    { 
                      reaction_types[i]=reaction_data.tau_type[reaction_index][i];
                      reaction_energies[i]=part_energy*reaction_data.tau_energy[reaction_index][i];
                      cout<<"tau decay with type is "<<reaction_types[i]<<" and reaction energy is "<<reaction_energies[i]<<endl;
                    }
                  }
                
                }
                part_energy*=reaction_data.tau_energy[reaction_index][0];
                //cout<<part_energy<<endl;
                //event.id[npart]=0;
                part_type=2;
              }
              if(part_type==4)
              {
                if(config.conversion)
                {
                for(int i=1;i<6;i++)
                  {
                    if(reaction_data.mu_type[reaction_index][i]!=-1) 
                    {
                      reaction_types[i]=reaction_data.mu_type[reaction_index][i];
                      reaction_energies[i]=part_energy*reaction_data.mu_energy[reaction_index][i];
                      cout<<"muon decay with type is "<<reaction_types[i]<<" and reaction energy is "<<reaction_energies[i]<<endl;
                    }
                  }
                
                }
                part_energy*=reaction_data.mu_energy[reaction_index][0];
                //event.id[npart]=2;
                part_type=1; 
              }

              // count the tau decay and the total number of interactions
              //event.ndk++;
              //event.npart++;
              //npart++;
            
              // Save the new particle energy and distance in the event structure/
              //event.E1[npart]=Energy_GeV;
              //event.v1[npart]=Ldist;
              
              // New particle is tagged as a tau neutrino.
              //event.id[npart]=0;
              //tag=0;
              // Comment next line to INCLUDE regeneration due to nu_tau(CC)->tau(decay)->nu_tau
              // break; // do not regenerate i.e. if nu_tau produced
            }
          }
            
          // Energy of tau below threshold => stop
          //if(Energy_GeV < Elim) break;
              
          // Flag particle below threshold
          //  brkcnt=1;
              
          // Count the number of times the particle fell below threshold
          //  brkcnt2++;
              
          //  break;}
            
          //} // end of 'if (tag == 1||3)' (particle is tau)

          //time to pop all the extra particles back in the stack to be looped over after
          if(config.conversion)
          {
            for( int i=1;i<6;i++)
            {
                
              if((reaction_energies[i]>Elim)&&(part_pos<Lmax2)&&(reaction_types[i]!=-1)) //particles inside earth and above threshold are stacked
              {
                cout<<"pushing new particle to stack \n";
                particle_data.part_type.push(reaction_types[i]);
                particle_data.part_energy.push(reaction_energies[i]);
                particle_data.part_pos.push(part_pos);
                particle_data.anti.push(anti);
                particle_data.generation.push(generation);
                particle_data.NC_num.push(NC_num);
                particle_data.CC_num.push(CC_num);
                particle_data.dc_num.push(dc_num);
                }
            }
          }
          if(part_energy<Elim)break; //uncomment to try to save partoic;e that fa;; below threshold
        } // ends loop 'for(Ldist=0.;Ldist<Lmax;)'
        // If the propagation was not terminated, then save the final energy of the particle in the event structure.
        //if(brkcnt!=1) event.Eend = Energy_GeV;
        
        //========================================
        // Tau above threshold emerging from Earth
        //========================================
        /*
        if((tag==1)&&(Energy_GeV>=Elim))
          {
          // Record this in the event structure.
          event.trig=1;
          
          // Generate a random number to sample shower position
          rndm = ((double) rand() / (double)(RAND_MAX));
          
          // Estimate the shower position (not used in this code)
          event.Shheight=(-1./dPdesdx(Energy_GeV,tag)*log(1.-rndm)+1000000.)*sin(PI*(refTheta/180.-1./2));
          event.Shlong=(-1./dPdesdx(Energy_GeV,tag)*log(1.-rndm)+1000000.)*cos(PI*(refTheta/180.-1./2));
          }
        */
        
        //=================================================
        // Write energy of emerging tau to output text file
        //=================================================
        //cout<<"save- type: "<<part_type<<" .energy: "<<part_energy<<".pos: "<<part_pos<<". gen: "<< generation<<endl;
      
        for(int i=0; i<4;i++)
        {
          //cout<<type_to_save[i]<<","<<part_type<<endl;
          //cout<<"part energy is "<<part_energy<<"with threshold "<<Elim<<endl;
          if((part_type==type_to_save[i])&&(part_energy>Elim))//changge to just poarticle tyoe if saving all particls even below thrwhpold
            {
              //cout<< "saving one particle of type "<<part_type<<endl;
              if (!config.energy_distribution) 
              {
                //cout<<"particlesaved!!!"<<endl;
                outEnergies << part_type<<" "<<NC_num << " " << CC_num << " " <<
                dc_num << " " << generation << " " <<
                setprecision(7)  << log10(part_energy)+9 << " " << log10(Energy_GeV)+9<<" "<<log10(Elim)+9<<endl;
              } 
              else 
              {
                outEnergies << part_type<<" "<< NC_num << " " << CC_num << " " <<
                dc_num << " " << generation << " " <<
                setprecision(7)  << part_energy << " " <<  //log10(event.Estart)+9.0 << " " << endl;
                endl;
              }
            }
          
        }//if the main particle isn't saved then it is forgetten falling below thrshold
      } while(!particle_data.part_type.empty()); // end of loop of stack
    } // end of loop for initial particles
      
    outEnergies << "END" << endl; // write END in the last line of the text output file.
    outEnergies.flush();
  
  }// end of loop ever all angles
  cout << "END" << endl;  // write END in the command line  
  
  return 0;
  
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// ===================================================
// Several functions
// ===================================================

//loads config.txt file
void load_config()
{
  ifstream fin("config.txt");
  string line;
  while (getline(fin,line))
  {
     istringstream sin(line.substr(line.find("=")+1));
     if (line.find("starting_type")!=-1) sin>>config.starting_type;
     else if (line.find("sim_angles")!=-1) sin>>config.sim_angles;
     else if (line.find("low_angle")!=-1) sin>>config.low_angle;
     else if (line.find("high_angle")!=-1) sin>>config.high_angle;
     else if (line.find("regen")!=-1) sin>>config.regen;
     else if (line.find("conversion")!=-1) sin>>config.conversion;
     else if (line.find("energy_distribution")!=-1) sin>>config.energy_distribution;
  }
}
// ########################################################
// CC neutrino cross-section (cm2) - various models fitted
// ########################################################

double dsigCC(double E, int CCmode, int type,int AntiNu )
{
    if(type==1||type==4)
    {
        double f=0.;
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

        // Connolly+, 2011 lower model (ARW's parametrization)
        double p0[4] = {-4.26355014e+01,   4.89151126e-01,   2.94975025e-02,  -1.32969832e-03};
        // Connolly+, 2011 middle model (ARW's parametrization)
        double p1[4] = { -5.35400180e+01,   2.65901551e+00, -1.14017685e-01,   1.82495442e-03};
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
    }
    if(type==2||type==5)
    {
        double f=0.;
        
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
        
        if(CCmode==0){
            // Connolly+, 2011 middle model (ARW's parametrization)
            double p[4] = { -5.35400180e+01,   2.65901551e+00, -1.14017685e-01,   1.82495442e-03};
            double log10_E_eV = log10(E)+9.; // E is in GeV, this parameterization is in log_10 ( E / eV ).
            for (int ii = 0 ; ii<4; ii++){
            f += p[ii]*pow(log10_E_eV,ii);
            //printf("\t %1.2e %1.2f %1.2e %1.2e %1.2e \n", E, log10_E_eV, p[ii], pow(log10_E_eV,ii), f);
            }
            f = pow(10,f);
            // printf("CC middle %1.2e %1.2f %1.2e %1.2f\n", E, log10_E_eV, f, log10(f));
        }
        
        if(CCmode==1){
            // Connolly+, 2011 lower model (ARW's parametrization)
            double p[4] = {-4.26355014e+01,   4.89151126e-01,   2.94975025e-02,  -1.32969832e-03};
            double log10_E_eV = log10(E)+9.; // E is in GeV, this parameterization is in log_10 ( E / eV ).
            for (int ii = 0 ; ii<4; ii++){
            f += p[ii]*pow(log10_E_eV,ii);
            //printf("\t %1.2e %1.2f %1.2e %1.2e %1.2e \n", E, log10_E_eV, p[ii], pow(log10_E_eV,ii), f);
            }
            f = pow(10,f);
            // printf("CC lower %1.2e %1.2f %1.2e %1.2f\n", E, log10_E_eV, f, log10(f));
        }
        
        if(CCmode==2){
            // Connolly+, 2011 upper model (ARW's parametrization)
            double p[4] = {-5.31078363e+01,   2.72995742e+00,  -1.28808188e-01,   2.36800261e-03};
            double log10_E_eV = log10(E)+9.; // E is in GeV, this parameterization is in log_10 ( E / eV ).
            for (int ii = 0 ; ii<4; ii++){
            f += p[ii]*pow(log10_E_eV,ii);
            //printf("\t %1.2e %1.2f %1.2e %1.2e %1.2e \n", E, log10_E_eV, p[ii], pow(log10_E_eV,ii), f);
            }
            f = pow(10,f);
            // printf("CC upper %1.2e %1.2f %1.2e %1.2f\n", E, log10_E_eV, f, log10(f));
        }
        
        return f;
    }
}

double dsigNC(double E, int CCmode, int type,int AntiNu)
{
    if(type==1||type==4)
    {
        double f=0.;
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

        // Connolly+, 2011 lower model (ARW's parametrization)
        double p0[4] = {-4.42377028e+01, 7.07758518e-01, 1.55925146e-02, -1.02484763e-03};
        // Connolly+, 2011 middle model (ARW's parametrization)
        double p1[4] = { -5.41463399e+01,   2.65465169e+00,  -1.11848922e-01,   1.75469643e-03};
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
    if(type==2||type==5)
    {
        double f=0.;
        
        // 	double l1=log10(E);
        // 	double l2=l1*l1;
        // 	double l3=l2*l1;
        // 	double l4=l3*l1;
        // 	double l5=l4*l1;
        // 	double l6=l5*l1;
        // 	double l7=l6*l1;
        // 	double l8=l7*l1;
        // 	double l9=l8*l1;
        // 	f = pNC0+pNC1*l1+pNC2*l2+pNC3*l3+pNC4*l4+pNC5*l5+pNC6*l6+pNC7*l7+pNC8*l8+pNC9*l9;
        
        //      f = 5.00969*pow(E,0.34944)*1e-36;
        
        // See references of models in dsigCC
        //      f = (-36.3345965603+7.14693605311*pow(E,0.293313250614))/2.4*1.e-36; // CKMT
        //	f = (-280.544665122+10.3452620208*pow(E,0.317119535055))/2.4*1.e-36; // ALLM
        //      f = (-799.252409182+52.4932827684*pow(E,0.244551044541))/2.4*1.e-36; // ASW
        // f = (-259.30822396+9.31732621406*pow(E,0.302056103343))*1.e-36;      // Sarkar
        
        if(CCmode==0){
            // Connolly+, 2011 middle model (ARW's parametrization)
            double p[4] = { -5.41463399e+01,   2.65465169e+00,  -1.11848922e-01,   1.75469643e-03};
            double log10_E_eV = log10(E)+9.; // E is in GeV, this parameterization is in log_10 ( E / eV ).
            for (int ii = 0 ; ii<4; ii++){
            f += p[ii]*pow(log10_E_eV,ii);
            //printf("\t %1.2e %1.2f %1.2e %1.2e %1.2e \n", E, log10_E_eV, p[ii], pow(log10_E_eV,ii), f);
            }
            f = pow(10,f);
            // printf("NC middle %1.2e %1.2f %1.2e %1.2f\n", E, log10_E_eV, f, log10(f));
        }
        
        if(CCmode==1){
            // Connolly+, 2011 lower model (ARW's parametrization)
            double p[4] = {-4.42377028e+01, 7.07758518e-01, 1.55925146e-02, -1.02484763e-03};
            double log10_E_eV = log10(E)+9.; // E is in GeV, this parameterization is in log_10 ( E / eV ).
            for (int ii = 0 ; ii<4; ii++){
            f += p[ii]*pow(log10_E_eV,ii);
            //printf("\t %1.2e %1.2f %1.2e %1.2e %1.2e \n", E, log10_E_eV, p[ii], pow(log10_E_eV,ii), f);
            }
            f = pow(10,f);
            // printf("NC lower %1.2e %1.2f %1.2e %1.2f\n", E, log10_E_eV, f, log10(f));
        }
        
        if(CCmode==2){
            // Connolly+, 2011 upper model (ARW's parametrization)
            double p[4] = {-5.36713302e+01,   2.72528813e+00,  -1.27067769e-01,   2.31235293e-03};
            double log10_E_eV = log10(E)+9.; // E is in GeV, this parameterization is in log_10 ( E / eV ).
            for (int ii = 0 ; ii<4; ii++){
            f += p[ii]*pow(log10_E_eV,ii);
            //printf("\t %1.2e %1.2f %1.2e %1.2e %1.2e \n", E, log10_E_eV, p[ii], pow(log10_E_eV,ii), f);
            }
            f = pow(10,f);
            // printf("NC upper %1.2e %1.2f %1.2e %1.2f\n", E, log10_E_eV, f, log10(f));
        }
        return f;
    }
}

// ###################################################
// 1./(gamma*c*muon0) with muon0 lifetime of muon (cm^-1)
// ###################################################
double dPdesdx(double E, int type)
{
    double f;
    if(type==1||type==4) f=mmuon/(E*muondl);
    if(type==2||type==5) f=mtau/(E*taudl);
        
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
  if(type==1||type==4)
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
  if(type==2||type==5)
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
  
  return f;
}

// ###################################################
// ###################################################
double funcalph(double *x, int *par, int type)
// double tfuncalph(double *x)
{
  if(type==2||type==5)
  {
    double f;
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
    return f;
  }
  if(type==1||type==4)
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
    return f;
  }

  
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
    if(type==2) //tau neutrino
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
    if(type==1) //muon neutrino
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
    string list_of[10000];
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
    
    for(int count=0;count<10000;count++)
    {
        getline(tau_decays,line);
        list_of[count]=line;
       
    }
    
    //parse string for products and fractional energies.
    for(int i=0;i<10000;i++)
    {
        for(int z=0;z<6;z++) {products[z]="-1,0.0";};
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
    for(int i=0;i<10000;i++)
    {
        for(int z=0;z<6;z++) {products[z]="-1,0.0";}
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

int convert_types(int pythia_type)
{
    int type=-1;
    switch (pythia_type){

        case -15:
            type=5;
            break;
        case 15:
            type=5;
            break;
        case -16:
            type=2;
            break;
        case 16:
            type=2;
            break;
        case 13:
            type=4;
            break;
        case -13:
            type=4;
            break;
        case 14:
            type=1;
            break;
        case -14:
            type=1;
            break;
        case 12:
            type=0;
            break;
        case -12:
            type=0;
            break;
        case 11:
            type=3;
            break;
        case -11:
            type=3;
            break;    
        default:
            type=-1;   
    }
    return type;
}



