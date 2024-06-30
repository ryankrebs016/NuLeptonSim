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
  //full sim things
  string data_dir;            //directory to store the data files
  int anti;                   //+1 is normal, -1 if anti particle
  int starting_type;          // starting type of neutrinos
  bool regen;                 // true or false for lepton regeneration through decays
  bool conversion;            // true or false to include produced particles from W +/- decays
  bool energy_distribution;   // true or false to use energy distribution
  bool use_sto_inst_of_cont;
  string min_muon_sto_loss;
  string min_tau_sto_loss;
  int n_throws;
  int n_traj;
  bool save_emerging;
  double energy_threshold;    //Set energy threshold

  //det things
  int detector;               // decides detector type, 0=Earth Emerging, 1=Spherical Volume, 2=Cylindrical Volume
  double det_volume;          //Det volume used to calculte radius for spherical model
  double ice_det_depth;
  double ice_det_rad;

  //save neutrino things
  bool save_neutrinos;        //Choose is neutrinos are saved
  bool save_nue;
  bool save_numu;
  bool save_nutau;

  //save charged lepton things
  bool save_charged;          //Choose is charged leptons are saved
  bool save_mu;
  bool save_tau;

  //save events in volume things
  bool save_events;           //choose to save events
  bool save_nu_events;
  bool save_sto_events;
  bool save_dec_events;

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

double decay_length(double P, double E, int particle_type);
bool in_volume(double x,double y,double z,double det_rad, double det_depth);
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
void generate_trajectory(double* xi, double* xf);

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