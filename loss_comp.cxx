#include"include/sto_losses.h"
#include"include/cont_losses.h"
//#include"include/Constants.h"
#include<stdio.h>
#include<math.h>

#include<iostream>
#include<fstream>
#include<ostream> 

#include<time.h>

using namespace std;

double dPdesdx(double E, int type)
{
    double f;
    if(type==13||type==14) f=mmuon/(E*muondl);
    if(type==15||type==16) f=mtau/(E*taudl);
    
    return f;
}

double L_decay(double P, double E, int particle_type)
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
    double decay_length=-log(1-P)*E*c_light*lifetime/mass;
    return decay_length;
}


void usage()
{
    cout<<"1st arg is initial energy\n second arg is type\n third arg is number to sim\n";
}


int main(int argc, char * argv[])
{

    srand(time (NULL));
    //ofstream out_rand("out_rand.txt");
    //for(int i=0;i<100000;i++)
    //{
    //    double rando=(double)rand()/(double)RAND_MAX*1000;
    //    out_rand<<rando<<endl;
    //}
    //out_rand.flush();
    //out_rand.close();
    //exit(1);


    //input args


    double i_energy=1E6; //units in GeV. so this i 1PeV
    cout<<argc<<endl;
    if(argc>1) 
    {
        i_energy=atof(argv[1]);
        cout<<"using energy "<<i_energy<<endl;
    }
    int type=15;
    if(argc>2)
    {
        type=atoi(argv[2]);
        cout<<"using type "<<type<<endl;
    }

    int num_parts=1000;
    if(argc>3)
    {
        num_parts=atof(argv[3]);
        cout<<"using num parts "<<num_parts<<endl;
    }

    printf("%i particles at 10^%.1f GeV of type %i\n",num_parts,log10(i_energy),type);

    continuous_loss_prop cont;
    stochastic_lepton_prop sto;
    sto.load_tables();


    double log_i_energy=log10(i_energy);
    double thresh=1E2; //GeV, lowest energy on LUT is 10^11 eV


   
    double dens=0.92; //ice
    //double dens=2.65; //rock

    sto.set_table_pointer(dens);
    time_t cont1;
    time_t cont2;
    time_t sto1;
    time_t sto2;
    cont.set_values(i_energy,dens,type,0);
    sto.set_val(i_energy,type,dens);
    
    int num_part=1;
    double cont_dist=0;
    double sto_dist=0;
    double energy=i_energy;
    int cont_loop=0;
    int sto_loop=0;
    int sto_decays=0;
    int cont_decays=0;

    double cont_dists[num_parts];
    double sto_dists[num_parts];
    double rando;double dL;
    cout<<"cont"<<endl;
    bool decays;
    double sto_e[num_parts];
    double cont_e[num_parts];
    double decay_length=0;

    //ofstream sto_loss("test_sto.txt");
    string data_dir="data/";

    time(&cont1);
    for(int i =0;i<num_parts;i++)
    {
        
        cont_dists[i]=0;
        cont_e[i]=0;
        if (i%(num_parts/10) ==0) cout<<i<<endl;
        cont.set_values(i_energy,dens,type,0);
        energy=i_energy;
        int dist_loop=0;
        while(energy>thresh)
        {
            dist_loop++;
            
            cont_loop++;
            dL=cont.get_interaction_length();
            cont_dists[i]+=dL;
            
            rando=(double)rand()/(double)RAND_MAX;
            
            //cout<< rando<<endl;
            //cout<<dL<<","<<dPdesdx(energy,type)<<endl;
            decays= rando<(1.-exp(-dL*dPdesdx(energy,type)));
            //cout<<dL<<endl;
            energy=energy-cont.get_energy_loss();
            //cout<<d<<endl;
            if (decays)
            {
                //cout<<"decay"<<endl;
                cont_decays++;
                
                break;
            }

            cont.set_values(energy,dens,type,0);

            
        
        }
        cont_dist+=cont_dists[i]/num_parts;
        cont_e[i]=energy;
        //cont.set_zero();
        
    }
    time(&cont2);
    energy=i_energy;
    time(&sto1);
    cout<<"sto"<<endl;
    for(int i=0;i<num_parts;i++)
    {
        //dens=0.92;
        sto_dists[i]=0;
        sto_e[i]=0;
        if (i%(num_parts/10) ==0) cout<<i<<endl;
        sto.set_val(i_energy,type,dens);
        energy=i_energy;
        while(energy>thresh)
        {
            sto_loop++;
            dL=sto.get_interaction_length();
            if(dL<0.)
            {
                cout<<"negative step, oospie"<<endl;
                exit(1);

            }
            
            double samp=sto.get_sampled_energy();
            //while(samp>.99807) 
            //{
                
            //    dL=sto.get_interaction_length();
            //    samp=sto.get_sampled_energy();
                //cout<<"how many times can I get stuck here"<<endl;

            //}
            
            //if(samp>.58)samp=.58; //this worked surprisingly well at matching dist. to cont
            //.6-> mean cont = 70.87km, mean sto = 68. 10^12GeV tau w/ 10^6 GeV threshold
            //matches really well at low energies
            //the high loss events are needed! to stay consistant to the cont. but the magnitude of them should
            //be lower to match... arb clamping at .6, 1/2 samp, might need to adjust cdf or the method the energy
            //get produced?

            //.58 best fit at 10^21 eV

            //samp=samp*2./3.;
            //sto_loss<<sto.sto_type<<","<<samp<<endl;
            decay_length=L_decay((double)rand()/(double)RAND_MAX,energy,type)*1e2; //m to cm
            //cout<<decay_length<<" vs "<<dL<<endl;
            if(decay_length<0)
            {
                cout<<"decay oopsie"<<endl;
                exit(-1);
            }
            rando=(double)rand()/(double)RAND_MAX;
            //decays= rando<(1.-exp(-dL*dPdesdx(energy,type)));

            if(decay_length>dL)
            {
                //no decay, do interaction step
                decays=false;
                sto_dists[i]+=dL;
                //cout<<dL<<endl;
                energy=(1-samp)*energy;
            }
            else //(decay_length<=dL)
            {
                decays=true;
                //do decay, no interaction step, do decay step
                sto_dists[i]+=decay_length;
            }

            //sto_dists[i]+=dL;
            //cout<<dL<<endl;
            //energy=(1-samp)*energy;

            //cout<<dL<<endl;
            //if(sto_dists[i]<1E6)dens=.92;
            //else if (sto_dists[i]>=1E6)dens=10.;
            

            
            //cout<<sto.get_sampled_energy()<<endl;
            if (decays)
            {
                sto_decays++;
                
                //cout<<"decays"<<endl;
                break;
            }

            sto.set_val(energy,type,dens);
        
        }
        sto_e[i]=energy;
        sto_dist+=sto_dists[i]/num_parts;
    }
    
    time(&sto2);
    
    //printf("cont dist %f km\n",cont_dist/pow(10,5));
    printf("sto dist %f km\n",sto_dist/100000);
    printf("cont dist %f km\n",cont_dist/100000);

    printf("cont looped %f time. sto looped %f time.\n",log10(cont_loop),log10(sto_loop));
    printf("cont took %li s and sto took %li s\n",(cont2-cont1),(sto2-sto1));
    printf("cont decays %i. sto decays %i\n",cont_decays,sto_decays);

    string cont_en_filename="data/cont_en";
    string cont_dist_filename="data/cont_dist";
    string sto_en_filename="data/sto_en";
    string sto_dist_filename="data/sto_dist";
    int padding=3;
    if (log_i_energy>9.99) padding=4;
    cont_en_filename+="_";
    cont_en_filename+=to_string(log_i_energy).substr(0,padding);
    cont_en_filename+="_";
    cont_en_filename+=to_string(type);
    cont_en_filename+=".txt";

    cont_dist_filename+="_";
    cont_dist_filename+=to_string(log_i_energy).substr(0,padding);
    cont_dist_filename+="_";
    cont_dist_filename+=to_string(type);
    cont_dist_filename+=".txt";

    sto_en_filename+="_";
    sto_en_filename+=to_string(log_i_energy).substr(0,padding);
    sto_en_filename+="_";
    sto_en_filename+=to_string(type);
    sto_en_filename+=".txt";

    sto_dist_filename+="_";
    sto_dist_filename+=to_string(log_i_energy).substr(0,padding);
    sto_dist_filename+="_";
    sto_dist_filename+=to_string(type);
    sto_dist_filename+=".txt";

    ofstream out_cont_en(cont_en_filename.c_str());
    ofstream out_sto_en(sto_en_filename.c_str());
    //sto_loss.flush();
    //sto_loss.close();


    ofstream out_cont(cont_dist_filename.c_str());
    ofstream out_sto(sto_dist_filename.c_str());
    
    cout<<"saving\n";
    for(int i=0;i<num_parts;i++)
    {
        if(i%(num_parts/10)==0)cout<<i<<endl;
        out_cont<<cont_dists[i]<<endl;
        out_sto<<sto_dists[i]<<endl;
        out_cont_en<<cont_e[i]<<endl;
        out_sto_en<<sto_e[i]<<endl;
    }

    out_sto_en.flush();
    out_sto_en.close();
    out_cont_en.flush();
    out_cont_en.close();
    out_cont.flush();
    out_cont.close();

    out_sto.flush();
    out_sto.close();
    
    return 0;

}


