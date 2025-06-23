
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
#include "Simu_elost.h"

using namespace std;
// ===================================================
// Several functions
// ===================================================

int still_going_towards_det(double * pos, double * traj_vec, double det_rad, double det_depth)
{
  //check if in volume
  bool in_vol = in_volume(pos[0], pos[1], pos[2], det_rad,det_depth);

  if (in_vol) return 1;
  else
  {
    //check if above top and if going up
    double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    //double elevation_angle = tan2(pos[2],r);

    //above going up
    if(pos[2] > R0 && traj_vec[2] > 0) return 2;

    //below going down
    if(pos[2] < R0-det_depth*1e5 && traj_vec[2] < 0) return 3;

    //outside going away from center
    double r_dot_traj = pos[0]*traj_vec[0]+pos[1]*traj_vec[1];//+pos[2]*traj_vec[2];
    if(r_dot_traj > 1e-9) return 4;

    //made it here is is still heading to the det
    return 0;
  }
}

void generate_trajectory(double* xi, double* xf, double rad, double depth, double exit_angle_cutoff)
{
  //xi starting point on the earth
  //xf ending point on the earth

  //vertex location in volume and cart directions
  double xv, yv, zv; 
  double dx, dy, dz;

  //get rand r, z, and ang
  double r = sqrt((double)rand()/(double)RAND_MAX)*rad*1e5;
  double z = (double)rand()/(double)RAND_MAX*depth*1e5;
  double ang = (double)rand()/(double)RAND_MAX*360*PI/180;

  //transform to cartesian
  xv = r*sin(ang);
  yv = r*cos(ang);
  zv = R0-z;

  //check if marginally above the sphere
  if(xv*xv+yv*yv+zv*zv > R02) 
  {
    //set to top of sphere
    zv = sqrt(R02-xv*xv-yv*yv)-10;//-100; 
  }

  double phi = (double)rand()/(double)RAND_MAX*360*PI/180;

  //convert exit angle to an elevation angle and take cos
  double cos_cutoff = sin((exit_angle_cutoff-90)*PI/180);
  double theta = (double)rand()/(double)RAND_MAX*(cos_cutoff+1)-1;
  theta = asin(theta);

  dx = cos(theta)*cos(phi);
  dy = cos(theta)*sin(phi);
  dz = sin(theta);

  //find intersection of neutrino line segment to earth surface
  double t_start;
  double t_end;

  double a = (dx*dx+dy*dy+dz*dz);
  double b = 2*(dx*xv+dy*yv+dz*zv);
  double c = xv*xv+yv*yv+zv*zv-R02;

  t_start = (-b - sqrt(b*b-4*a*c))/(2*a);
  t_end = (-b + sqrt(b*b-4*a*c))/(2*a);

  //set entrance and exit points
  double buff = .999999; //really depends how much padding is needed, ex 1 vs .9999999
  xi[0] = (xv+dx*t_start*buff);
  xi[1] = (yv+dy*t_start*buff);
  xi[2] = (zv+dz*t_start*buff);

  xf[0] = (xv+dx*t_end*buff);
  xf[1] = (yv+dy*t_end*buff);
  xf[2] = (zv+dz*t_end*buff);
  
}

int get_intersections(double *pos, double * dir, double radius, double * int_dists)
{
  //math math
  double d_start;
  double d_end;

  double a = (dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
  double b = 2*(dir[0]*pos[0]+dir[1]*pos[1]+dir[2]*pos[2]);
  double c = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]-radius*radius;
  double desc = b*b-4*a*c;

  if (desc >= 0)
  {
    d_start = (-b - sqrt(desc))/(2*a);
    d_end = (-b + sqrt(desc))/(2*a);

    //two solutions
    int_dists[0] = d_start;
    int_dists[1] = d_end;
    
    if (desc < 1e-9) return 1;
    else return 2;
  }
  else
  {
    //no solutions
    int_dists[0]=0;
    int_dists[1]=0;
    return 0;
  }
}
double analytic_distance_from_grammage(double *pos, double *dir, double grammage, double new_depth, double new_dens)
{
  //initialize layers (for now 8 layers but could be increased to n)
  int num_layers = 10;
  double layer_radii[10] = {122150000,348000000,570100000,577100000,597100000,615100000,634660000,635600000,637800000,637800001}; //cm
  double layer_mean_density[10] = {12.884,10.898,4.904,3.958,3.837,3.490,3.351,2.9,2.6,new_dens};

  //assign new layer vals
  if(new_depth) layer_radii[8] = 637800000-new_depth*1e5;
  if(new_dens) layer_mean_density[9] = new_dens;

  //intialize vars
  double temp_pos[3] = {pos[0], pos[1], pos[2]};
  double next_pos[3] = {0.};
  double traversed_length = 0;
  double radius = sqrt(temp_pos[0]*temp_pos[0]+temp_pos[1]*temp_pos[1]+temp_pos[2]*temp_pos[2]);
  double current_grammage = 0;
  float dir_dot_r = 0;
  bool going_up = false;
  bool left_earth = false;
  int current_layer = 9;
  double layer_grammage = 0;
  double lengths[2] = {0.};
  int num_sols = 0;
  double min_bump = 1; //cm... really just a safety thing
  int search_layer = 0;
  int which_distance = 0;
  
  //search inside out since the depths aren't nicely spaced
  for (int i=0; i<10; i++)
  {
    if(radius < layer_radii[i])
    {
      current_layer = i;
      break;
    }
  }

  //loop until grammage distance found or left earth (should only happen at most 2*n-1 times)
  while(current_grammage<grammage || left_earth)
  {
    radius = sqrt(temp_pos[0]*temp_pos[0]+temp_pos[1]*temp_pos[1]+temp_pos[2]*temp_pos[2]);
    
    //check if dir is currently pointed towards center of earth or away
    dir_dot_r = -temp_pos[0]*dir[0]-temp_pos[1]*dir[1]-temp_pos[2]*dir[2];
    if (dir_dot_r<0) going_up = true; 
    else going_up = false;

    if(going_up)
    {
      //in layer and points up so going to intersect its own surface only
      //get path length and add grammage
      num_sols = get_intersections(temp_pos, dir, layer_radii[current_layer], lengths);
      layer_grammage = layer_mean_density[current_layer]*lengths[1];

      if(current_grammage+layer_grammage > grammage)
      {
        //hits grammage before leaving
        //maybe update position if doing inside func
        return traversed_length+(grammage-current_grammage)/layer_mean_density[current_layer]; 
      }

      else
      {
        //leaves earth
        current_grammage += layer_grammage;
        traversed_length += layer_grammage/layer_mean_density[current_layer];
        lengths[1] += min_bump;
        temp_pos[0] = temp_pos[0]+dir[0]*lengths[1];
        temp_pos[1] = temp_pos[1]+dir[1]*lengths[1];
        temp_pos[2] = temp_pos[2]+dir[2]*lengths[1];

        if(current_layer == num_layers-1)
        {
          left_earth = true;
          current_layer += 1;
          return traversed_length+1e9;
        }
        current_layer += 1;
        continue;
      }
    }

    else
    {
      // path points down so maybe lower layer or maybe skimming
      // determine if in core and adjust indexes
      if(current_layer == 0) 
      {
        search_layer = current_layer;
        which_distance = 1;
      }
      else 
      {
        search_layer = current_layer-1;
        which_distance = 0;
      }

      num_sols = get_intersections(temp_pos, dir, layer_radii[search_layer], lengths);

      if(num_sols > 0)
      {
        // travels to lower layer or exits core
        // going down so going to have two valid solutions ahead of it. take "behind" as first point
        layer_grammage = layer_mean_density[current_layer]*lengths[which_distance];

        if(current_grammage+layer_grammage > grammage)
        {
          // hits grammage target
          // maybe update position if doing inside func
          return traversed_length+(grammage-current_grammage)/layer_mean_density[current_layer]; 
        }
        else
        {
          // goes down a layer
          current_grammage += layer_grammage;
          traversed_length += layer_grammage/layer_mean_density[current_layer];
          lengths[which_distance] += min_bump;
  
          temp_pos[0] = temp_pos[0]+dir[0]*lengths[which_distance];
          temp_pos[1] = temp_pos[1]+dir[1]*lengths[which_distance];
          temp_pos[2] = temp_pos[2]+dir[2]*lengths[which_distance];

          current_layer -= 1;
          continue;
        }
      }

      else
      {
        //skims through layer and goes up a layer or exit earth
        //passes through layer without going down so use current layer radius
        num_sols = get_intersections(temp_pos, dir, layer_radii[current_layer], lengths);

        //get distance to current layer intersection
        layer_grammage = layer_mean_density[current_layer]*lengths[1];

        if(current_grammage+layer_grammage > grammage)
        {
          //hits grammage before leaving
          //maybe update position if doing inside func
          return traversed_length+(grammage-current_grammage)/layer_mean_density[current_layer]; 
        }
        else
        {
          //leaves earth
          current_grammage += layer_grammage;
          traversed_length += layer_grammage/layer_mean_density[current_layer];
          lengths[1] += min_bump;
  
          temp_pos[0] = temp_pos[0]+dir[0]*lengths[1];
          temp_pos[1] = temp_pos[1]+dir[1]*lengths[1];
          temp_pos[2] = temp_pos[2]+dir[2]*lengths[1];

          if(current_layer == num_layers-1)
          {
            left_earth = true;
            return traversed_length+1e9;
          }
          current_layer += 1;
          continue;
        }
      }
    }
  }
}

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
    printf("wrong particle\n");
    exit(-1);
  }

  if(1-P<0)
  {
    printf("log bug\n");
    exit(-1);
  } 

  double decay_length = -log(1-P)*E*c_light*lifetime/mass*1e2;
  return decay_length; //cm
}

void set_points_from_angle(input_file *input, double angle)
{
  double x_step = cos((angle-90)*PI/180);
  double y_step = sin((angle-90)*PI/180);
  
  double start_point[2];
  double end_point[2];
  double depth = 0;
  double small_r = 0;
  end_point[0] = 0;
  end_point[1] = R0-depth+small_r*y_step;
  printf("End point %f, %f\n",end_point[0],end_point[1]);

  double point_slope=  0;
  double y_pos = 0;
  double y_neg = 0;

  if(angle!=180.0 || angle!=0)
  {
    point_slope = -tan((angle-90)*PI/180);
    double det = 4*(R0-depth)*(R0-depth)-4*(1+point_slope*point_slope)*((R0-depth)*(R0-depth)-R02*point_slope*point_slope);
    double denom = 2*(1+point_slope*point_slope);
    if(det<0)
    {
      printf("bad determinant\n");
      //return -1;
    }
    y_pos = (2*(R0-depth)+sqrt(det))/denom;
    y_neg = (2*(R0-depth)-sqrt(det))/denom;
  }

  if(angle >= 90.0) start_point[1] = y_neg;
  if(angle < 90.0) start_point[1] = y_pos;
  
  start_point[0] = sqrt(R02-start_point[1]*start_point[1]);
  printf("Start point %f, %f\n",start_point[0],start_point[1]);

  if(angle==180.0)
  {
    start_point[1] = -R0;
    start_point[0] = 0;
  }
  if(angle==0.)
  {
    start_point[1] = R0;
    start_point[0] = 0;
  }
  input[0].eex = 0;
  input[0].eey = 0;
  input[0].eez = R0;

  input[0].xe = start_point[0];
  input[0].ye = 0;
  input[0].ze = start_point[1];

  input[0].xi = start_point[0];
  input[0].yi = 0;
  input[0].zi = start_point[1];

  input[0].xf = 0;
  input[0].yf = 0;
  input[0].zf = R0;
}

bool in_volume(double x, double y, double z, double det_rad, double det_depth)
{
  //cylinder
  double rad=sqrt(x*x+y*y);
  if(rad<=det_rad*1e5 && (z<R0) && (z>R0-det_depth*1e5))
  {
    return true;
  }
  else return false;

  /* sphere
  if(detec.inner_rad==0)
  {
   double sp_R=sqrt(x*x+y*y+z*z);
   if(sp_R<detec.inner_sphere) return true;
   else return false;
  }
  */
}

int load_input(input_file *in, int file_length, string filename)
{
  fstream test_csv(filename.c_str());
  string row = "";
  getline(test_csv,row);
  
  double holding[12];
  int row_ind = 0;
  double conv = pow(10,5);//km to cm
  while(row.size()!=0)
  {
    int i = 0;
    int place = 0;
    int length = 0;
    int index = 0;
    while(row.find(",",place)!=-1)
    {
      
      index = row.find(",", place);
      length = index-place;
      holding[i] = conv*atof(row.substr(place,length).c_str());
      place = index+1;
      i++;

    }

    holding[11] = conv*atof(row.substr(place,row.size()-place).c_str());
    
    in[row_ind].xi = holding[0];
    in[row_ind].yi = holding[1];
    in[row_ind].zi = holding[2];
    in[row_ind].eex = holding[3];
    in[row_ind].eey = holding[4];
    in[row_ind].eez = holding[5];
    in[row_ind].xf = holding[9];
    in[row_ind].yf = holding[10];
    in[row_ind].zf = holding[11];
    in[row_ind].xe = holding[6];
    in[row_ind].ye = holding[7];
    in[row_ind].ze = holding[8];

    getline(test_csv, row);

    row_ind++;
    if(row_ind+1 > file_length) 
    {
      printf("Reached eof, %i > %i\n",row_ind,file_length);
      return 0;
    }
  }
  return 0;
}

string make_particle_dir(int argc, char **argv, string out_dir, string es_temp, string angs_temp, config_t * config)
{
  string nameEnergies = "";
  nameEnergies += config->data_dir;
  nameEnergies += "/particles";
  nameEnergies += "/particles_";
  nameEnergies += es_temp;
  nameEnergies += "_";
  nameEnergies += angs_temp;
  
  if(argc==8)
  {
    nameEnergies += "_";
    nameEnergies += argv[6];
    nameEnergies += "km_ice_";
    if(atof(argv[4])==0) nameEnergies += "mid";
    if(atof(argv[4])==1) nameEnergies += "low";
    if(atof(argv[4])==2) nameEnergies += "upp";
    nameEnergies += "CS_";
    if(atof(argv[5])==0) nameEnergies += "std";
    if(atof(argv[5])==1) nameEnergies +=" low";
    nameEnergies += "EL";
  };
  nameEnergies+=".dat";
  
  return nameEnergies;
}
string make_event_dir(int argc, char **argv, string out_dir, string es_temp, string angs_temp, int p_type, string label, config_t * config)
{
  string type_temp = "";
  switch (p_type)
  {
    case 12:
      type_temp="nue";
      break;
    case 14:
      type_temp="numu";
      break;
    case 16:
      type_temp="nutau";
      break;
    case -12:
      type_temp="anti_nue";
      break;
    case -14:
      type_temp="anti_numu";
      break;
    case -16:
      type_temp="anti_nutau";
      break;
    case 0:
      type_temp="mixed";
      break;
  }
  string nameEvents = "";
  nameEvents += config->data_dir;
  nameEvents += "/events/";
  nameEvents += type_temp;
  nameEvents += "_events_";
  nameEvents += es_temp;
  if (label!="")
  {
    nameEvents += "_";
    nameEvents += label;
  }

  //nameEvents+="_";
  //nameEvents+=angs_temp;
  nameEvents += ".dat";

  return nameEvents;

}

string make_lepton_dir(int argc, char **argv, string out_dir, string es_temp, string angs_temp, int p_type, config_t * config)
{
  string type_temp="";
  switch (p_type)
  {
    case 12:
      type_temp="nue";
      break;
    case 14:
      type_temp="numu";
      break;
    case 16:
      type_temp="nutau";
      break;
    case -12:
      type_temp="anti_nue";
      break;
    case -14:
      type_temp="anti_numu";
      break;
    case -16:
      type_temp="anti_nutau";
      break;
    case 0:
      type_temp = "mixed";
      break;
  }

  string nameLeptons = "";
  nameLeptons += config->data_dir;
  nameLeptons += "/leptons/";
  nameLeptons += type_temp;
  nameLeptons += "_leptons_";
  nameLeptons += es_temp;
  nameLeptons += "_";
  nameLeptons += angs_temp;
  nameLeptons += ".dat";
  return nameLeptons;

}

double get_dens_from_coords(double *coords, Earth * terra)
{
  double f;
  double radius = sqrt(coords[0]*coords[0]+coords[1]*coords[1]+coords[2]*coords[2]);
  f = terra->GetDensity(radius);
  return f;
}

//loads config.txt file
void load_config(config_t * config)
{
  ifstream fin("config.txt");
  string line;
  while (getline(fin,line))
  {
    if(line.find("#")!=-1) continue;
     istringstream sin(line.substr(line.find("=")+1));
     if((int)line.find("data_dir")!=-1) sin>>config->data_dir;
     else if ((int)line.find("starting_type")!=-1) sin>>config->starting_type;
     else if ((int)line.find("anti")!=-1) sin >>config->anti;
     else if ((int)line.find("regen")!=-1) sin>>config->regen;
     else if ((int)line.find("conversion")!=-1) sin>>config->conversion;
     else if ((int)line.find("energy_distribution")!=-1) sin>>config->energy_distribution;
     else if ((int)line.find("detector")!=-1) sin>>config->detector;
     else if ((int)line.find("det_volume")!=-1) sin >>config->det_volume;
     else if ((int)line.find("save_neutrinos")!=-1) sin >>config->save_neutrinos;
     else if ((int)line.find("save_charged")!=-1) sin >>config->save_charged;
     else if ((int)line.find("energy_threshold")!=-1) sin >>config->energy_threshold;
     else if ((int)line.find("save_events")!=-1) sin >>config->save_events;
     else if ((int)line.find("run_throws")!=-1) sin >>config->run_throws;
     else if ((int)line.find("n_throws")!=-1) sin >>config->n_throws;
     else if ((int)line.find("run_number")!=-1) sin >>config->run_number;
     else if ((int)line.find("num_emerging_leptons")!=-1) sin >>config->num_emerging_leptons;
     else if ((int)line.find("n_traj")!=-1) sin >>config->n_traj;
     else if ((int)line.find("save_nu_events")!=-1) sin>>config->save_nu_events;
     else if ((int)line.find("save_sto_events")!=-1) sin>>config->save_sto_events;
     else if ((int)line.find("save_dec_events")!=-1) sin>>config->save_dec_events;
     else if ((int)line.find("save_emerging")!=-1) sin>>config->save_emerging;
     else if ((int)line.find("use_sto_inst_of_cont")!=-1) sin>>config->use_sto_inst_of_cont;
     else if ((int)line.find("min_muon_sto_loss")!=-1) sin>>config->min_muon_sto_loss;
     else if ((int)line.find("min_tau_sto_loss")!=-1) sin>>config->min_tau_sto_loss;     
     else if ((int)line.find("default_energy")!=-1) sin>>config->default_energy;
     else if ((int)line.find("default_angle")!=-1) sin>>config->default_angle;
     else if ((int)line.find("default_cc")!=-1) sin>>config->default_cc;
     else if ((int)line.find("default_eloss")!=-1) sin>>config->default_eloss;
     else if ((int)line.find("default_layer_thickness")!=-1) sin>>config->default_layer_thickness;
     else if ((int)line.find("default_layer_density")!=-1) sin>>config->default_layer_density;
     else if ((int)line.find("ice_det_depth")!=-1) sin>>config->ice_det_depth;
     else if ((int)line.find("ice_det_rad")!=-1) sin>>config->ice_det_rad;
     else if ((int)line.find("save_nue")!=-1) sin>>config->save_nue;
     else if ((int)line.find("save_numu")!=-1) sin>>config->save_numu;
     else if ((int)line.find("save_nutau")!=-1) sin>>config->save_nutau;
     else if ((int)line.find("save_mu")!=-1) sin>>config->save_mu;
     else if ((int)line.find("save_tau")!=-1) sin>>config->save_tau;
     else if ((int)line.find("sto_force_distance")!=-1) sin>>config->sto_force_distance;
     else if ((int)line.find("ext_traj")!=-1) sin>>config->ext_traj;
     else if ((int)line.find("ang_cutoff")!=-1) sin>>config->ang_cutoff;
     else if ((int)line.find("save_final_part_state")!=-1) sin>>config->save_final_part_state;
  }
}

void load_geo(det_t * det) //could be useful later on
{
  ifstream fin("det_geom.txt");
  string line;
  while (getline(fin,line))
  {
     istringstream sin(line.substr(line.find("=")+1));
     if ((int)line.find("inner_rad")!=-1) sin>>det->inner_rad;
     else if ((int)line.find("inner_depth")!=-1) sin >>det->inner_depth;
     else if ((int)line.find("outer_rad")!=-1) sin>>det->outer_rad;
     else if ((int)line.find("outer_depth")!=-1) sin>>det->outer_depth;
     else if ((int)line.find("outer_rad")!=-1) sin>>det->outer_rad;
     else if ((int)line.find("outer_depth")!=-1) sin>>det->outer_depth;
  }
  if(det->outer_rad==0) det->traj_weights = (2*det->outer_rad*PI*PI + det->outer_depth*det->outer_rad*2*PI)/R02;
  if(det->outer_sphere==0) det->traj_weights = 2*PI*det->outer_sphere*det->outer_sphere/R02;
  
}
// ########################################################
// CC neutrino cross-section (cm2) - various models fitted
// ########################################################

double dsigGR(double E, int type, int AntiNu)
{
  //check if it is anti nue
  if(type!=12 || AntiNu!=-1) return 0.;

  double GF = 1.166E-5; //GeV^-2
  double A = 3.4986E-27;
  double gammaW = 2.085; //GeV/c^2

  double cross_section_num = A*4*GF*GF*E*mW*mW*mW*mW*me;
  double cross_section_denom = 3*2*PI*((mW*mW-2*me*E)*(mW*mW-2*me*E)+mW*mW*gammaW*gammaW);

  if(cross_section_denom==0) printf("GR cross section div 0 err.\n");

  double cs = cross_section_num/cross_section_denom;

  return cs;
}

double dsigCC(double E, int CCmode, int type,int AntiNu )
{
  
  double f=0.;
  double p[4];
  
  AntiNu = int(abs(AntiNu-1)/2);  //return anti (+1 as part., -1 as anti part.) to 1 = antiparticle
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

  if( E < E_switch )
  {
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

    
  // Unused parameterizations

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


double dsigNC(double E, int CCmode, int type, int AntiNu)
{
  double f=0.; 
  AntiNu = int(abs(AntiNu-1)/2); //return anti (+1 as part., -1 as anti part.) to 1 = antiparticle
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

  if( E < E_switch )
  {
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

    if(type==13 || type==14) 
      f = mmuon/(E*muondl);
    else //(type==15 || type==16) 
      f = mtau/(E*taudl);

    return f;
}


// #######################################################
// Mean Earth density along a chord at angle theta
// #######################################################
double mean_dens_chord(double theta, Earth * terra)
{
  double f;
  double z = 0.;
  f = mean(&theta, &z, terra);
  return f;
}


//// ###################################################
//// Mean Earth density along a chord of length Lmax
//// ###################################################
////mac double mean(double *x, double *par)
double mean(double *x, double *par, Earth * terra)
{
  double Lmax = 2*R0*cos(PI-PI*x[0]/180.);
  
  // stupid numerical integral
  double sum = 0.;
  int N = 1000000;
  double dx = Lmax/((double) (N));
  for (int ii=0; ii<=N; ii++)
  {
    double x_val = ((double) ii)*dx;
    sum += earthdens(&x_val, &Lmax, terra) * dx;
  }
  sum /= Lmax;

  return sum;
}

// #######################################################
// Earth density as a function of radius - read from table
// #######################################################
double earthdens( double *x, double *par, Earth * terra)
{
  double f;
  double Current_Radius = sqrt(R02 - (x[0]*par[0])+x[0]*x[0]);
  f = terra->GetDensity(Current_Radius);

  return f;
}				

//initializes the reaction data arrays from the pythia tables
void initialize_reaction(int tau_type[][6], int mu_type[][6], double tau_ene[][6], double mu_ene[][6])
{
  
    string line;
    string list_of[100000];
    string reaction;
    string particle[6];
    string energy[6];
    string delimiter_p = ":";
    string delimiter_e = ",";
    string products[6];
    int counter;
    int done;

    //taus
    ifstream tau_decays("tables/pythia_tau.txt");
    
    for(int count=0; count<100000; count++)
    {
        getline(tau_decays,line);
        list_of[count] = line;
       
    }
    
    //parse string for products and fractional energies.
    for(int i=0; i<100000; i++)
    {
      for(int z=0; z<6; z++) {products[z] = "0,0.0";}
      counter = 0;
      done = 0;
        
      for(int j=0; j<6; j++)
      {
        if(list_of[i].find(delimiter_p) != string::npos)
        {
          products[j] = list_of[i].substr(0,list_of[i].find(delimiter_p));
          list_of[i].erase(0, list_of[i].find(delimiter_p)+1);
          counter++;
        }
        if(list_of[i].find(delimiter_p) == string::npos && done==0)
        {
          products[counter] = list_of[i];
          done=1;
        }
      }
        
      for(int s=0; s<6; s++)
      {
        tau_type[i][s] = convert_types(stoi(products[s].substr(0,products[s].find(delimiter_e))));
        products[s].erase(0, products[s].find(delimiter_e)+1);
        tau_ene[i][s] = stod(products[s]);
      }
    }
    tau_decays.close();

    //muons
    ifstream mu_decays("tables/pythia_muon.txt");
    
    for(int count=0; count<10000; count++)
    {
      getline(mu_decays,line);
      list_of[count]=line;
    }

    //parse string for products and fractional energies.
    for(int i=0; i<100000; i++)
    {
      for(int z=0; z<6; z++) {products[z] = "0,0.0";}
      counter = 0;
      done = 0;
      for(int j=0; j<6; j++)
      {
        if(list_of[i].find(delimiter_p) != string::npos)
        {
          products[j] = list_of[i].substr(0,list_of[i].find(delimiter_p));
          list_of[i].erase(0, list_of[i].find(delimiter_p)+1);
          counter++;
        }
        if(list_of[i].find(delimiter_p) == string::npos && done==0)
        {
          products[counter] = list_of[i];
          done = 1;
        }
      }
        
      for(int s=0; s<6; s++)
      {
        mu_type[i][s] = convert_types(stoi(products[s].substr(0,products[s].find(delimiter_e))));
        products[s].erase(0, products[s].find(delimiter_e)+1);
        mu_ene[i][s] = stod(products[s]);
      }
    }
  mu_decays.close();
}

void make_dirs(string dirs)
{
  mkdir(dirs.c_str(), 0755);
  mkdir((dirs+"/particles").c_str(), 0755);
  mkdir((dirs+"/events").c_str(), 0755);
  mkdir((dirs+"/LUT").c_str(), 0755);
  mkdir((dirs+"/leptons").c_str(), 0755);
}

int convert_types(int pythia_type) 
{
  int type=0;
  switch (pythia_type)
  {
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

  return type;
}