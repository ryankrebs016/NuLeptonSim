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
  double default_energy;
  double default_angle;
  int default_cc;
  int default_eloss;
  double default_layer_thickness;
  double default_layer_density;
  bool run_throws;
  int n_throws;
  bool run_number;
  int num_emerging_leptons;
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
  bool save_final_part_state;
  bool ext_traj;
  string trajectory_file;
  double ang_cutoff;
  double sto_force_distance;
}config_init;   // data struct to hold the value read from config file



typedef struct
{
  //input file parameters setting startpoint, end point, and etrance point to the volume
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
}input_file;


typedef struct
{
  //cylindrical
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


int still_going_towards_det(double * pos, double * traj_vec, double det_rad, double det_depth);
bool in_volume(double x,double y,double z,double det_rad, double det_depth);

double decay_length(double P, double E, int particle_type);

string make_particle_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp);
string make_event_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp, int p_type, string label);
string make_lepton_dir(int argc, char **argv,string out_dir,string es_temp,string angs_temp, int p_type);

// initialize reaction from pythia table and convert pythia type to tag code used in code
void initialize_reaction(int tau_type[][6], int mu_type[][6],double tau_ene[][6],double mu_ene[][6]);
int convert_types(int pythia_type);

//load values from config file
void load_config();
void load_geo();
int load_input(input_file *in,int length, string filename);

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
// Local density as a function of x,y,z coordinates
double get_dens_from_coords(double *coords);

// -------------------------------------------------
// Average density as a function of zenith angle
double mean( double *x, double *par);
double mean_dens_chord(double theta);

//---------------------------------------------------------------
void make_dirs(string dirs);
void generate_trajectory(double* xi, double* xf,double rad, double depth, double exit_ang_cutoff);
void set_points_from_angle(input_file *input,double angle);
double analytic_distance_from_grammage(double *pos, double *dir, double grammage, double new_depth, double new_dens);
int get_intersections(double *pos, double * dir, double radius, double * int_dists);
