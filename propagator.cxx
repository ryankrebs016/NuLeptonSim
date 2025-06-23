
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
#include <vector>

#include "Table.hh"
#include "Earth.hh"
#include "Constantes.hh"

#include "sto_losses.h"
#include "cont_losses.h"
#include "Simu_elost.h"

//#define LUT

using namespace std;

reaction_tables_t reaction_data;
config_t config;
det_t det;
Earth *terra = new Earth(0.0, 2.6); 

typedef struct
{
    int id;
    double pos[3];
    double dir[3];
    double energy;
    int type;
    int anti;
} part_info_t;

typedef struct
{
    int id;
    double pos[3];
    double dir[3];
    double particle_energy;
    double shower_energy;
    int particle_type;
    int interaction_type;

} event_info_t;

int input_translator(string line, part_info_t * part_info)
{
    //place holder values for now
    part_info->id = 0;
    part_info->type = 13;
    part_info->pos[0] = 0;
    part_info->pos[1] = -30*1e5;
    part_info->pos[2] = 6378*1e5-2*1e5;
    part_info->dir[0] = 0;
    part_info->dir[1] = 1;
    part_info->dir[2] = 0;
    part_info->energy = 1e10; //Gev
    part_info->anti=1;

    return 0;
}
string output_translator(event_info_t * evt)
{
  printf("Event: id %i, type %i, energy %.2f ",evt->id,evt->particle_type, evt->particle_energy);
  printf("x %.2f, y %.2f, z %.2f ",evt->pos[0],evt->pos[1], evt->pos[2]);
  printf("dx %.2f, dy %.2f, dz %.2f ",evt->dir[0],evt->dir[1], evt->dir[2]);
  printf("shower type %i, shower energy %.2f\n",evt->interaction_type,evt->shower_energy);
  
  return "";
}

int print_part(part_info_t * part)
{
  printf("Particle: id %i, type %i, energy %.2f ",part->id,part->type*part->anti, part->energy);
  printf("x %.2f, y %.2f, z %.2f ",part->pos[0],part->pos[1], part->pos[2]);
  printf("dx %.2f, dy %.2f, dz %.2f\n",part->dir[0],part->dir[1], part->dir[2]);
  return 0;
}

int main(int argc, char **argv)
{
  double time_start = time(NULL); //set start time

  // load config and charged lepton energy loss classes
  load_config(&config); 
  stochastic_lepton_prop sto;
  sto.load_tables(config.min_muon_sto_loss, config.min_tau_sto_loss);
  
  for(int i=0; i<100000; i++)
  {
    for(int j=0; j<6; j++)
    {
      reaction_data.tau_type[i][j] = reaction_data.mu_type[i][j] = 0;
      reaction_data.tau_energy[i][j] = reaction_data.mu_energy[i][j] = 0.0;
    }
  }  
  initialize_reaction(reaction_data.tau_type,reaction_data.mu_type,reaction_data.tau_energy,reaction_data.mu_energy);
  
  terra->depth_new_layer = config.default_layer_thickness;
  terra->dens_new_layer = config.default_layer_density;

  //string input_filename = argv[1];
  //string output_filename = argv[2];

  double threshold_energy = config.energy_threshold * 1e-9;

  // open input and output file

  string line;
  string outline;
  int count = 0;
  //while(getline(input_file,line))
  printf("det: rad %f km, depth %f km\n",config.ice_det_rad,config.ice_det_depth);
  while(count<10000)
  {
    //if(count% int(.1*count)) printf("count %i\n",count);
    count++;
    double dens = 0;
    double dL = 0;
    part_info_t part_info = {0.};
    event_info_t event_info = {0.};
    input_translator(line, &part_info);
    
    int going_towards_det = still_going_towards_det(part_info.pos,part_info.dir,config.ice_det_rad, config.ice_det_depth);
    double sphere_radius2 = part_info.pos[0]*part_info.pos[0]+part_info.pos[1]*part_info.pos[1]+part_info.pos[2]*part_info.pos[2];
    dens=get_dens_from_coords(part_info.pos, terra);
    double traversed_grammage=0;
    bool break_loop = false;
    print_part(&part_info);
    //printf("%i\n",going_towards_det);
    while(going_towards_det <2 && sphere_radius2 < R02 && !break_loop)
    {
      
      //print_part(&part_info);
      if(abs(part_info.type) != 13 && abs(part_info.type) != 15)
      dens=get_dens_from_coords(part_info.pos, terra);
      double frac_loss = 0;
      double sampled_decay_length = 0;
      double random_num = 0;

      int reaction_types[6] = {0,0,0,0,0,0};
      double reaction_energies[6] = {0,0,0,0,0,0};
      int anti_type[6] = {1,1,1,1,1,1};

      //get a decay length
      random_num = (double)rand()/(double)RAND_MAX;
      sampled_decay_length = decay_length((double)rand()/(double)RAND_MAX,part_info.energy,part_info.type); //cm

      // get the interaction info
      // force stochastic if in-ice det and below force distance set in config, or only use stochastic losses
      
      sto.set_val(part_info.energy, part_info.type, dens);
      dL = sto.get_interaction_length();
      frac_loss = sto.get_sampled_energy();

        //continuous
        //cont.set_values(part_energy, dens, part_type, 0);
        //dL = cont.get_interaction_length();
        //frac_loss = cont.get_energy_loss()/part_energy;
      

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
      
      if(sampled_decay_length > dL)
      {
        //==============================
        // The tau lepton does NOT decay
        //=============================
        // cout<<"no decay"<<endl;

        // check to save the deposition
        if(config.save_sto_events && in_volume(part_info.pos[0],part_info.pos[1],part_info.pos[2],config.ice_det_rad,config.ice_det_depth) && ((frac_loss*part_info.energy)>threshold_energy))
        {
          int temp_show = 0;
          if(sto.sto_type == 0) temp_show=1;
          if(sto.sto_type == 1) temp_show=1;//sto_type==0 brem, sto_type==1 pp, sto_type==2 pn
          if(sto.sto_type == 2) temp_show=0; 

          //outline = output_translater(event_info)
          //outputfile.write(outline);
          event_info.id = part_info.id;
          event_info.pos[0]=part_info.pos[0];
          event_info.pos[1]=part_info.pos[1];
          event_info.pos[2]=part_info.pos[2];
          event_info.dir[0]=part_info.dir[0];
          event_info.dir[1]=part_info.dir[1];
          event_info.dir[2]=part_info.dir[2];
          event_info.particle_energy = part_info.energy;
          event_info.shower_energy = part_info.energy * frac_loss;
          event_info.interaction_type = temp_show;
          event_info.particle_type = part_info.type * part_info.anti;
          output_translator(&event_info);
          //outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<<Energy_GeV << ","<<part_energy<<","<<frac_loss<<","
          //  <<part_type*anti<<","<<sto.sto_type+4<<","<<temp_show<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<sto_num<<"\n";
          //outLep<<log10(part_energy)<<","<<frac_loss<<","<<log10(part_energy*frac_loss)<<","<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<sto.sto_type<<","<<i<<","<<num_count<<","<<part_type<<","<<sto_num<<endl;
          //sto_num++;
        }

        // update position and energy
        // part_pos = part_pos+dL;
        part_info.pos[0] = part_info.pos[0]+dL*part_info.dir[0];
        part_info.pos[1] = part_info.pos[1]+dL*part_info.dir[1];
        part_info.pos[2] = part_info.pos[2]+dL*part_info.dir[2];
        part_info.energy = (1.-frac_loss)*part_info.energy;

        // Calculate the traversed grammage
        traversed_grammage+=dL*dens;
      }
      else
      {
        //=======================
        // The lepton decays
        //=======================
        // tau or muon does decay

        //if(in_volume(pos[0],pos[1],pos[2],config.ice_det_rad,config.ice_det_depth))
        //{
        //  decay_num++;
        //}

        //part_pos = part_pos+sampled_decay_length;
        part_info.pos[0] = part_info.pos[0]+sampled_decay_length*part_info.dir[0];
        part_info.pos[1] = part_info.pos[1]+sampled_decay_length*part_info.dir[1];
        part_info.pos[2] = part_info.pos[2]+sampled_decay_length*part_info.dir[2];


        // Calculate the traversed grammage
        traversed_grammage+=sampled_decay_length*dens;
      
        //dc_num++;
        //cout<<"tau path length before decaying is "<<tau_path<<endl;
        //cout<<"tau energy at decay "<<part_energy<<endl;
        //cout<<"dL is "<<dL<<endl;
        
        // Get the energy of the neutrino produced in the decay
        //generation++;

        //if(tag==1) TauData.ThrowFinal(finalstate);
        //if(tag==3) MuonData.ThrowFinal(finalstate);
        //Energy_GeV=finalstate[0]*Energy_GeV;

        double initial_energy = part_info.energy;
        int reaction_index = (double)rand()/(double)RAND_MAX*100000;

        int initial_particle_type = part_info.type;
        double dec_inel = 0;
        int initial_particle = part_info.type;
        int is_had_or_em = -1;
        if(part_info.type == 15)
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
                anti_type[j] = part_info.anti*anti_ness;
                reaction_types[j] = abs(reaction_data.tau_type[reaction_index][j]);
                reaction_energies[j] = part_info.energy*reaction_data.tau_energy[reaction_index][j];
              }
              // get what kind of shower will be produced !!!double check
              if(is_had_or_em == -1 && reaction_types[j] == 0) is_had_or_em = 0;//for HAD
              if(is_had_or_em == -1 && reaction_types[j] == 11) is_had_or_em = 1;//for EM
              if(is_had_or_em == -1 && reaction_types[j] == 13) is_had_or_em = 2;//for HAD+EM
            }
          }
          part_info.energy *= reaction_data.tau_energy[reaction_index][0];
          //part_info.type = 16;
          //part_info.anti = -part_info.anti;
        }
        if(part_info.type == 13)
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
                anti_type[j] = part_info.anti*anti_ness;
                reaction_types[j] = abs(reaction_data.mu_type[reaction_index][j]);
                reaction_energies[j] = part_info.energy*reaction_data.mu_energy[reaction_index][j];
              }
              if(is_had_or_em==-1 && reaction_types[j]==0) is_had_or_em=0;//for HAD
              if(is_had_or_em==-1 && reaction_types[j]==11) is_had_or_em=1;//for EM
            }
          }
          part_info.energy *= reaction_data.mu_energy[reaction_index][0];
          //part_info.type = 14;
          //part_info.anti = -part_info.anti; 
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
        if((initial_energy*frac_energy_dumped>threshold_energy) && config.save_dec_events && in_volume(part_info.pos[0],part_info.pos[1],part_info.pos[2],config.ice_det_rad,config.ice_det_depth))
        {

          event_info.id = part_info.id;
          event_info.pos[0]=part_info.pos[0];
          event_info.pos[1]=part_info.pos[1];
          event_info.pos[2]=part_info.pos[2];
          event_info.dir[0]=part_info.dir[0];
          event_info.dir[1]=part_info.dir[1];
          event_info.dir[2]=part_info.dir[2];
          event_info.particle_energy = part_info.energy;
          event_info.shower_energy = part_info.energy * frac_loss;
          event_info.interaction_type = is_had_or_em;
          event_info.particle_type = part_info.type * part_info.anti;

          output_translator(&event_info);
          // save the decay event
          //outEvents<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<x_step<<","<<y_step<<","<<z_step<<","<<initial_flavor<<","<<Energy_GeV << ","<<initial_energy<<","<<frac_energy_dumped<<","
          //    <<initial_particle*anti<<","<<3<<","<<is_had_or_em<<","<<NC_num<<","<<dc_num<<","<<GR_num<<","<<i<<","<<num_count<<","<<-1<<"\n";
          //event_count++;
          //out_gram<<i<<","<<num_count<<","<<part_pos<<","<<traversed_grammage<<","<<1<<endl;
          //got_event=true;
        }

        // dont track neutrino heres
        break_loop = true;
      }
      going_towards_det = still_going_towards_det(part_info.pos,part_info.dir,config.ice_det_rad, config.ice_det_depth);
      sphere_radius2 = part_info.pos[0]*part_info.pos[0]+part_info.pos[1]*part_info.pos[1]+part_info.pos[2]*part_info.pos[2];

    }
  }
  return 0;
}