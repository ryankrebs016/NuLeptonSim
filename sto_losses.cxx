#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
//#include "include/Constantes.hh"
#include <string>
#include <math.h>
#include "math.h"
#include <cstring>
#include <iomanip>
#include<cmath>
#include<time.h>

using namespace std;
//class to handle all stochastic losses in NuLeptonSim. requires initialization before particle loops
//and requires particle properties to be updated before getting sampled energy and interaction length

class stochastic_lepton_prop
{
    public:
    double energy_GeV;
    double log_energy_MeV;
    int p_type;
    int loss_mode;
    int loss_model;
    int sto_type;
    double dens;
    double 

    static const int ts1=51;
    static const int ts2=10000;
    double log_energies_MeV[ts1];
   
  
    //muon tables
    double cs_brem_muon[ts1];
    double cs_pp_muon[ts1];
    double cs_pn_muon[ts1];
    double cdf_val_muon[ts2];
    double **cdf_xs_brem_muon=new double *[ts1];
    double **cdf_xs_pp_muon=new double *[ts1];
    double **cdf_xs_pn_muon=new double *[ts1];
    double total_cs_muon[ts1];



    //tau tables
    double cs_brem_tau[ts1];
    double cs_pp_tau[ts1];
    double cs_pn_tau[ts1];
    double cdf_val_tau[ts2];
    double **cdf_xs_brem_tau=new double* [ts1];
    double **cdf_xs_pp_tau= new double *[ts1];
    double **cdf_xs_pn_tau= new double *[ts1];
    double total_cs_tau[ts1];

    void set_val(double energy,int particle,double temp_dens);
    void set_model(int temp_loss_mode, int temp_loss_model);
    double get_interaction_length();//uses loss_mode and loss_model and energy and density
    double get_sampled_energy();
    void load_tables(); //load tables into memory
    double interpolate_cs(int sto);
    double interpolate_int_length();
    void save_event(ofstream * lepton_file,double dep_energy);
    void set_sto_type();
};

void stochastic_lepton_prop::load_tables()
{   int time1=time(NULL);
    cout<<"loading tables"<<endl;;
    for(int i = 0;i<ts1;i++)
    {
        cdf_xs_brem_tau[i]=new double [ts2];
        cdf_xs_pp_tau[i]=new double [ts2];
        cdf_xs_pn_tau[i]=new double [ts2];
        cdf_xs_brem_muon[i]=new double [ts2];
        cdf_xs_pp_muon[i]=new double [ts2];
        cdf_xs_pn_muon[i]=new double [ts2];
    }
    
    string temp_str,temp_str2,temp_str3,temp_str4,temp_str5,temp_str6;
    string delim=" ";
    ifstream energy_file("stochastic_tables_ice/MeV_energies.txt");

    ifstream muon_cs_brem("stochastic_tables_ice/cs_brem_muon.txt");
    ifstream muon_cs_pp("stochastic_tables_ice/cs_pp_muon.txt"); 
    ifstream muon_cs_pn("stochastic_tables_ice/cs_pn_muon.txt");
    ifstream muon_cdf_val("stochastic_tables_ice/cdf_values_muon.txt");
    ifstream muon_cdf_xs_brem("stochastic_tables_ice/cdf_xs_brem_muon.txt");
    ifstream muon_cdf_xs_pp("stochastic_tables_ice/cdf_xs_pp_muon.txt");
    ifstream muon_cdf_xs_pn("stochastic_tables_ice/cdf_xs_pn_muon.txt");

    ifstream tau_cs_brem("stochastic_tables_ice/cs_brem_tau.txt");
    ifstream tau_cs_pp("stochastic_tables_ice/cs_pp_tau.txt"); 
    ifstream tau_cs_pn("stochastic_tables_ice/cs_pn_tau.txt");
    ifstream tau_cdf_val("stochastic_tables_ice/cdf_values_tau.txt");
    ifstream tau_cdf_xs_brem("stochastic_tables_ice/cdf_xs_brem_tau.txt");
    ifstream tau_cdf_xs_pp("stochastic_tables_ice/cdf_xs_pp_tau.txt");
    ifstream tau_cdf_xs_pn("stochastic_tables_ice/cdf_xs_pn_tau.txt");

    //printf("files open\n");
    for(int i=0;i<ts1;i++)
    {
        //cout<<i<<",";
        getline(energy_file,temp_str);
        log_energies_MeV[i]=atof(temp_str.c_str());
        //cout<<MeV_energy[i]<<"-";
        getline(muon_cs_brem,temp_str);
        cs_brem_muon[i]=atof(temp_str.c_str());
        
        getline(muon_cs_pp,temp_str);
        cs_pp_muon[i]=atof(temp_str.c_str());
        
        getline(muon_cs_pn,temp_str);
        cs_pn_muon[i]=atof(temp_str.c_str());
        
        getline(tau_cs_brem,temp_str);
        cs_brem_tau[i]=atof(temp_str.c_str());
        
        getline(tau_cs_pp,temp_str);
        cs_pp_tau[i]=atof(temp_str.c_str());
        
        getline(tau_cs_pn,temp_str);
        cs_pn_tau[i]=atof(temp_str.c_str());

        total_cs_muon[i]=cs_brem_muon[i]+cs_pp_muon[i]+cs_pn_muon[i];
        total_cs_tau[i]=cs_brem_tau[i]+cs_pp_tau[i]+cs_pn_tau[i];
        
        int counter=0;

        if(i==0)
        {
            for(int q=0;q<ts2;q++)
            {
                getline(muon_cdf_val,temp_str);
                cdf_val_muon[q]=atof(temp_str.c_str());
                //cout<<cdf_val_muon[j]<<",";
                getline(tau_cdf_val,temp_str);
                cdf_val_tau[q]=atof(temp_str.c_str());
            }
            //cout<<j<<" ";
        }

        getline(muon_cdf_xs_brem,temp_str);
        getline(muon_cdf_xs_pp,temp_str2);
        getline(muon_cdf_xs_pn,temp_str3);

        getline(tau_cdf_xs_brem,temp_str4);
        getline(tau_cdf_xs_pp,temp_str5);
        getline(tau_cdf_xs_pn,temp_str6);

        for(int j=0;j<ts2;j++)
        {
          
            if(temp_str.find(delim)!=string::npos)
            {
                cdf_xs_brem_muon[i][j]=atof(temp_str.substr(0,temp_str.find(delim)).c_str());
                temp_str.erase(0,temp_str.find(delim)+1);

                cdf_xs_pp_muon[i][j]=atof(temp_str2.substr(0,temp_str2.find(delim)).c_str());
                temp_str2.erase(0,temp_str2.find(delim)+1);

                cdf_xs_pn_muon[i][j]=atof(temp_str3.substr(0,temp_str3.find(delim)).c_str());
                temp_str3.erase(0,temp_str3.find(delim)+1);

                cdf_xs_brem_tau[i][j]=atof(temp_str4.substr(0,temp_str4.find(delim)).c_str());
                temp_str4.erase(0,temp_str4.find(delim)+1);

                cdf_xs_pp_tau[i][j]=atof(temp_str5.substr(0,temp_str5.find(delim)).c_str());
                temp_str5.erase(0,temp_str5.find(delim)+1);

                cdf_xs_pn_tau[i][j]=atof(temp_str6.substr(0,temp_str6.find(delim)).c_str());
                temp_str6.erase(0,temp_str6.find(delim)+1);

                
            }
            if(temp_str.find(delim)==string::npos)
            {
                cdf_xs_brem_muon[i][ts2-1]=atof(temp_str.c_str());
                cdf_xs_pp_muon[i][ts2-1]=atof(temp_str2.c_str());
                cdf_xs_pn_muon[i][ts2-1]=atof(temp_str3.c_str());
                cdf_xs_brem_tau[i][ts2-1]=atof(temp_str4.c_str());
                cdf_xs_pp_tau[i][ts2-1]=atof(temp_str5.c_str());
                cdf_xs_pn_tau[i][ts2-1]=atof(temp_str6.c_str());
                
                
            }


        }
        //.same format just like 16 times rip lol
        //and one more loop in here for the 2d arrays
    }
    energy_file.close();
    muon_cs_brem.close();
    muon_cs_pp.close();
    muon_cdf_val.close();
    muon_cdf_xs_brem.close();
    muon_cdf_xs_pp.close();
    muon_cdf_xs_pn.close();

    tau_cs_brem.close();
    tau_cs_pp.close();
    tau_cdf_val.close();
    tau_cdf_xs_brem.close();
    tau_cdf_xs_pp.close();
    tau_cdf_xs_pn.close();

    

   
    //printf("\n\n\n");
    //for(int i=0;i<ts1;i++)cout<<MeV_energy[i]<<",";
    //cout<<endl;
    //for(int i=0;i<ts2;i++)cout<<cdf_xs_brem_muon[0][i]<<",";
    //cout<<setprecision(9);
    //for(int i=0;i<ts2;i++)cout<<cdf_xs_brem_muon[0][i]<<",";
    cout<<"loaded tables "<<time(NULL)-time1<<" s"<<endl;

};

double stochastic_lepton_prop::interpolate_cs(int sto)
{
    int index= (log_energy_MeV - 5. )/(10./(ts1-1));
    if(index==ts1-1)
    {
    if (p_type == 13)
    {
        switch (sto)
        {
            case 0: return cs_brem_muon[index];
            case 1: return cs_pp_muon[index];
            case 2: return cs_pn_muon[index];
        }
            
    }
    if (p_type == 15)
    {
        switch (sto)
        {
            case 0: return cs_brem_tau[index];
            case 1: return cs_pp_tau[index];
            case 2: return cs_pn_tau[index];
        }
            
    }
    }
    else
    {
    if (p_type == 13)
    {
        switch (sto)
        {
            case 0: return (log_energies_MeV[index+1]-log_energy_MeV)*(cs_brem_muon[index])/.2+(log_energy_MeV-log_energies_MeV[index])*cs_brem_muon[index+1]/.2;
            case 1: return (log_energies_MeV[index+1]-log_energy_MeV)*(cs_pp_muon[index])/.2+(log_energy_MeV-log_energies_MeV[index])*cs_pp_muon[index+1]/.2;
            case 2: return (log_energies_MeV[index+1]-log_energy_MeV)*(cs_pn_muon[index])/.2+(log_energy_MeV-log_energies_MeV[index])*cs_pn_muon[index+1]/.2;
            
            //case 0: return log_energy_MeV*(cs_brem_muon[index+1]-cs_brem_muon[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+cs_brem_muon[index];
            //case 1: return log_energy_MeV*(cs_pp_muon[index+1]-cs_pp_muon[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+cs_pp_muon[index];
            //case 2: return log_energy_MeV*(cs_pn_muon[index+1]-cs_pn_muon[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+cs_pn_muon[index];
           
        }
            
    }
    if (p_type == 15)
    {
        switch (sto)
        {
            case 0: return (log_energies_MeV[index+1]-log_energy_MeV)*(cs_brem_tau[index])/.2+(log_energy_MeV-log_energies_MeV[index])*cs_brem_tau[index+1]/.2;
            case 1: return (log_energies_MeV[index+1]-log_energy_MeV)*(cs_pp_tau[index])/.2+(log_energy_MeV-log_energies_MeV[index])*cs_pp_tau[index+1]/.2;
            case 2: return (log_energies_MeV[index+1]-log_energy_MeV)*(cs_pn_tau[index])/.2+(log_energy_MeV-log_energies_MeV[index])*cs_pn_tau[index+1]/.2;
            //case 0: return log_energy_MeV*(cs_brem_tau[index+1]-cs_brem_tau[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+cs_brem_tau[index];
            //case 1: return log_energy_MeV*(cs_pp_tau[index+1]-cs_pp_tau[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+cs_pp_tau[index];
            //case 2: return log_energy_MeV*(cs_pn_tau[index+1]-cs_pn_tau[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+cs_pn_tau[index];
           
        }
            
    }
    }
    //should never reach this location but if it does
    printf("failed to interpolate cs for leptons bc of type\n");
    exit(EXIT_FAILURE);
}

double stochastic_lepton_prop::interpolate_int_length()
{
    //cout<<log10(energy_eV)<<endl;
    int index= (log_energy_MeV - 5. )/(10./(ts1-1));
    //cout<<index<<endl;
    //cout<<index<<endl;
    //cout<<"got here"<<endl;
    if(index==ts1-1)
    {
        if(p_type==13)return 22.*1.66E-24/total_cs_muon[index];
        if(p_type==15)return 22.*1.66E-24/total_cs_tau[index];
    }
    else
    {
        //cout<<"in here"<<endl;
        //cout<<MeV_energy[index]<<endl;
        //cout<<total_cs_tau[index+1]<<","<<total_cs_tau[index]<<endl;
        //cout<<MeV_energy[index+1]<<","<<MeV_energy[index]<<endl;
        //cout<<log10(energy_eV*pow(10,-6))<<endl;
        //cout<<total_cs_tau[index+1]<<","<<total_cs_tau[index]<<endl;

        //if(p_type==13) return 22.*1.66E-24*((1./total_cs_muon[index+1]-1./total_cs_muon[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+1./total_cs_muon[index])*log_energy_MeV;
        //if(p_type==15) return 22.*1.66E-24*((1./total_cs_tau[index+1]-1./total_cs_tau[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])+1./total_cs_tau[index])*log_energy_MeV;

        double f1=0;
        double f2=0;
        double val=0;

        if(p_type==13) 
        {
            f1=1./total_cs_muon[index+1];
            f2=1./total_cs_muon[index];

            

            //return 22.*1.66E-24*((1./total_cs_muon[index+1]-1./total_cs_muon[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])*log_energy_MeV+1./total_cs_muon[index]);
        
        }
        
        if(p_type==15) 
        {
            f1=1./total_cs_tau[index+1];
            f2=1./total_cs_tau[index];
            
            
            //return 22.*1.66E-24*((1./total_cs_tau[index+1]-1./total_cs_tau[index])/(log_energies_MeV[index+1]-log_energies_MeV[index])*log_energy_MeV+1./total_cs_tau[index]);
        }
        val=22.*1.66E-24*((log_energies_MeV[index+1]-log_energy_MeV)/0.2*f2+(log_energy_MeV-log_energies_MeV[index])/0.2*f1);
     
        return val;
    }
    
    //again, shouldnt reach this, buttt
    printf("failed to interpolate total cross section bc of type\n");
    exit(EXIT_FAILURE);
}

void stochastic_lepton_prop::set_val(double energy, int particle, double temp_dens)
{
    energy_GeV=energy;
    log_energy_MeV=log10(energy_GeV)+3.;
    p_type=particle;
    dens=temp_dens;
    if((p_type!=15&&p_type!=13)||dens<0.)
    {
        cout<<"wrong particle type"<<endl;
        exit(1);
    }
}

void stochastic_lepton_prop::set_model(int temp_loss_mode, int temp_loss_model)
{
    //empty unless i converge continuous losses here
}

double stochastic_lepton_prop::get_interaction_length()//needs work
{  
    
    //double rando=(double)rand()/(double)RAND_MAX;
    //cout<<rando<<endl;
    double rand_lin=(double)rand()/(double)RAND_MAX;//*M_E;
    //cout<<rand_lin<<endl;
    //double rand_lin2=(double)rand()/(double)RAND_MAX;//*M_E;
    //cout<<rand_lin<<","<<rand_lin2<<endl;
    //cout<<rand_lin<<endl;
    //double random_exponential;
    //if((1-rand_lin)==0) random_exponential=
    
    double random_exponential=-log(1.-rand_lin); //need to come up with a way to generate random exponentials like numpy.random.exponentioal
    //cout<<random_exponential<<endl;
    //double val=interpolate_int_length();
    //cout<<val<<endl;
    //cout<<interpolate_int_length()<<endl;
    //cout<<random_exponential*interpolate_int_length()/dens<<endl;
    
    return random_exponential*interpolate_int_length()/dens;

}

void stochastic_lepton_prop::set_sto_type()
{
    double brem,pp,pn,tot;
    brem=interpolate_cs(0);
    pp=interpolate_cs(1);
    pn=interpolate_cs(2);
    
    tot=brem+pp+pn;
    double ran=(double)rand()/(double)RAND_MAX;
    //printf("brem cs, pp cs, pn cs - %f,%f,%f\n",brem/tot,pp/tot,pn/tot);
    if(ran<=(brem/tot)) sto_type=0;
    else if(ran<=((brem+pp)/tot)) sto_type=1;
    else if (ran<=1.) sto_type=2;
    //cout<<sto_type<<endl;

}
    
double stochastic_lepton_prop::get_sampled_energy()
{
    set_sto_type();
    //cout<<sto_type<<endl;
    //2d interpolation based on energy and on a random number between 0,ts2
    

    
    int index1= (log_energy_MeV - 5. )/(10./(ts1-1));
    double rand_index=(double)rand()/(double)RAND_MAX * (double)(ts2-1.);//old way in lin space
    rand_index=rand_index/(ts2-1.)*log10(ts2);
    
    //lets leave this mess for when I'm caffeinated

    rand_index=pow(10,rand_index)-1.;//new way in log space
    double rando=(double)rand()/(double)RAND_MAX;
    //cout<<rando<<",";
    long int temp=rando*1E10;//sets precision to 10 decimal places
    //cout<<temp<<",";
    rando=temp/1E10;
    //cout<<rando<<",";
    int index2=-1000.*log10(1.-rando);

    
    double f00,f01,f10,f11;
    double sampled_energy;
    double f_top,f_bot;
    double d_val;
    
    if(p_type==13)
    {
        switch(sto_type)
        {
            case 0:
            {
                if(index1==ts1-1 && index2==ts2-1)
                {
                    sampled_energy=cdf_xs_brem_muon[index1][index2];
                    return sampled_energy;
                }
                else if(index1==ts1-1 && index2!=ts2-1)
                {
                    d_val=(cdf_val_muon[index2+1]-cdf_val_muon[index2]);
                    
                    f00=cdf_xs_brem_muon[index1][index2];
                    f10=cdf_xs_brem_muon[index1][index2+1];
                    sampled_energy=(cdf_val_muon[index2+1]-rando)/d_val*f00+(rando-cdf_val_muon[index2])/d_val*f10;

                    //sampled_energy=(f10-f00)/ts2*rand_index+f00;
                    return sampled_energy;
                    //1d interpolation along rand sample axis
                }
                else if(index1!=ts1-1 && index2==ts2-1)
                {
                    f00=cdf_xs_brem_muon[index1][index2];
                    f01=cdf_xs_brem_muon[index1+1][index2];
                    //sampled_energy=(f01-f00)/10.*log_energy_MeV+f00;
                    sampled_energy=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    return sampled_energy;
                    //1d interpolated along energy axis
                }
                else if(index1!=ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_muon[index2+1]-cdf_val_muon[index2];
                    f00=cdf_xs_brem_muon[index1][index2];
                    f01=cdf_xs_brem_muon[index1+1][index2];
                    f10=cdf_xs_brem_muon[index1][index2+1];
                    f11=cdf_xs_brem_muon[index1+1][index2+1];
                    //interpolate_sampled_energy()
                    //full 2d interpolation 
                    
                    //f_bot=(f01-f00)/10. * log_energy_MeV+f00;
                    //f_top=(f11-f10)/10. * log_energy_MeV+f10;
                    //sampled_energy=(f_top-f_bot)/ts2*rand_index+f_bot;

                    f_bot=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    f_top=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f10+(log_energy_MeV-log_energies_MeV[index1])/.2*f11;
                    sampled_energy=(cdf_val_muon[index2+1] - rando)/d_val*f_bot+(rando-cdf_val_muon[index2])/d_val*f_top;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                }

            }
            case 1:
            {
                if(index1==ts1-1 && index2==ts2-1)
                {
                    sampled_energy=cdf_xs_pp_muon[index1][index2];
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                }
                else if(index1==ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_muon[index2+1]-cdf_val_muon[index2];
                    f00=cdf_xs_pp_muon[index1][index2];
                    f10=cdf_xs_pp_muon[index1][index2+1];
                    
                    sampled_energy=(cdf_val_muon[index2+1]-rando)/d_val*f00+(rando-cdf_val_muon[index2])/d_val*f10;
                    //sampled_energy=(index2+1-rand_index)/1.*f00+(rand_index-index2)/1.*f10;
                    //sampled_energy=(f10-f00)/ts2*rand_index+f00;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                    //1d interpolation along rand sample axis
                }
                else if(index1!=ts1-1 && index2==ts2-1)
                {
                    f00=cdf_xs_pp_muon[index1][index2];
                    f01=cdf_xs_pp_muon[index1+1][index2];
                    //sampled_energy=(f01-f00)/10.*log_energy_MeV+f00;
                    sampled_energy=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                    //1d interpolated along energy axis
                }
                else if(index1!=ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_muon[index2+1]-cdf_val_muon[index2];
                    f00=cdf_xs_pp_muon[index1][index2];
                    f01=cdf_xs_pp_muon[index1+1][index2];
                    f10=cdf_xs_pp_muon[index1][index2+1];
                    f11=cdf_xs_pp_muon[index1+1][index2+1];
                    //cout<<"f00,f01: "<<f00<<","<<f01<<endl;
                    //cout<<"f10,f11: "<<f10<<","<<f11<<endl;
                    //interpolate_sampled_energy()
                    //full 2d interpolation 
                    f_bot=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    f_top=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f10+(log_energy_MeV-log_energies_MeV[index1])/.2*f11;
                    //cout<<"f_bot,f_top: "<<f_bot<<","<<f_top<<endl;
                    //sampled_energy=(index2+1 - rand_index)/1.*f_bot+(rand_index-index2)/1.*f_top;
                    sampled_energy=(cdf_val_muon[index2+1] - rando)/d_val*f_bot+(rando-cdf_val_muon[index2])/d_val*f_top;

                    
                   
                

                   

                    //f_bot=(f01-f00)/10. * log_energy_MeV+f00;
                    //f_top=(f11-f10)/10. * log_energy_MeV+f10;
                    //sampled_energy=(f_top-f_bot)/ts2*rand_index+f_bot;
                    
                  
                    return sampled_energy;
                }
               
            }
            case 2:
            {
                
                if(index1==ts1-1 && index2==ts2-1)
                {
                    sampled_energy=cdf_xs_pn_muon[index1][index2];
                    return sampled_energy;
                }
                else if(index1==ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_muon[index2+1]-cdf_val_muon[index2];

                    f00=cdf_xs_pn_muon[index1][index2];
                    f10=cdf_xs_pn_muon[index1][index2+1];

                    //sampled_energy=(f10-f00)/ts2*rand_index+f00;
                    sampled_energy=(cdf_val_muon[index2+1]-rando)/d_val*f00+(rando-cdf_val_muon[index2])/d_val*f10;
                    //sampled_energy=(index2+1-rand_index)/1.*f00+(rand_index-index2)/1.*f10;
                    return sampled_energy;
                    //1d interpolation along rand sample axis
                }
                else if(index1!=ts1-1 && index2==ts2-1)
                {
                    f00=cdf_xs_pn_muon[index1][index2];
                    f01=cdf_xs_pn_muon[index1+1][index2];
                    //sampled_energy=(f01-f00)/10.*log_energy_MeV+f00;
                    sampled_energy=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    return sampled_energy;
                    //1d interpolated along energy axis
                }
                else if(index1!=ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_muon[index2+1]-cdf_val_muon[index2];
                    f00=cdf_xs_pn_muon[index1][index2];
                    f01=cdf_xs_pn_muon[index1+1][index2];
                    f10=cdf_xs_pn_muon[index1][index2+1];
                    f11=cdf_xs_pn_muon[index1+1][index2+1];
                    //interpolate_sampled_energy()
                    //full 2d interpolation 
                    
                    //f_bot=(f01-f00)/10. * log_energy_MeV+f00;
                    //f_top=(f11-f10)/10. * log_energy_MeV+f10;
                    //sampled_energy=(f_top-f_bot)/ts2*rand_index+f_bot;

                    f_bot=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    f_top=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f10+(log_energy_MeV-log_energies_MeV[index1])/.2*f11;
                    sampled_energy=(cdf_val_muon[index2+1] - rando)/d_val*f_bot+(rando-cdf_val_muon[index2])/d_val*f_top;

                    //sampled_energy=(index2+1 - rand_index)/1.*f_bot+(rand_index-index2)/1.*f_top;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                }
            }
        }
    }
    else if(p_type==15)
      {
        switch(sto_type)
        {
            case 0:
            {
                
                if(index1==ts1-1 && index2==ts2-1)
                {
                    sampled_energy=cdf_xs_brem_tau[index1][index2];
                    return sampled_energy;
                }
                else if(index1==ts1-1 && index2!=ts2-1)
                {
                    d_val=(cdf_val_tau[index2+1]-cdf_val_tau[index2]);
                    
                    f00=cdf_xs_brem_tau[index1][index2];
                    f10=cdf_xs_brem_tau[index1][index2+1];
                    sampled_energy=(cdf_val_tau[index2+1]-rando)/d_val*f00+(rando-cdf_val_tau[index2])/d_val*f10;

                    //sampled_energy=(f10-f00)/ts2*rand_index+f00;
                    return sampled_energy;
                    //1d interpolation along rand sample axis
                }
                else if(index1!=ts1-1 && index2==ts2-1)
                {
                    f00=cdf_xs_brem_tau[index1][index2];
                    f01=cdf_xs_brem_tau[index1+1][index2];
                    //sampled_energy=(f01-f00)/10.*log_energy_MeV+f00;
                    sampled_energy=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    return sampled_energy;
                    //1d interpolated along energy axis
                }
                else if(index1!=ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_tau[index2+1]-cdf_val_tau[index2];
                    f00=cdf_xs_brem_tau[index1][index2];
                    f01=cdf_xs_brem_tau[index1+1][index2];
                    f10=cdf_xs_brem_tau[index1][index2+1];
                    f11=cdf_xs_brem_tau[index1+1][index2+1];
                    //interpolate_sampled_energy()
                    //full 2d interpolation 
                    
                    //f_bot=(f01-f00)/10. * log_energy_MeV+f00;
                    //f_top=(f11-f10)/10. * log_energy_MeV+f10;
                    //sampled_energy=(f_top-f_bot)/ts2*rand_index+f_bot;

                    f_bot=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    f_top=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f10+(log_energy_MeV-log_energies_MeV[index1])/.2*f11;
                    sampled_energy=(cdf_val_tau[index2+1] - rando)/d_val*f_bot+(rando-cdf_val_tau[index2])/d_val*f_top;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                }

                //filler_cs=&cdf_xs_brem_tau[0][0];
            }
            case 1:
            {
                
                if(index1==ts1-1 && index2==ts2-1)
                {
                    sampled_energy=cdf_xs_pp_tau[index1][index2];
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                }
                else if(index1==ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_tau[index2+1]-cdf_val_tau[index2];
                    f00=cdf_xs_pp_tau[index1][index2];
                    f10=cdf_xs_pp_tau[index1][index2+1];
                    
                    sampled_energy=(cdf_val_tau[index2+1]-rando)/d_val*f00+(rando-cdf_val_tau[index2])/d_val*f10;
                    //sampled_energy=(index2+1-rand_index)/1.*f00+(rand_index-index2)/1.*f10;
                    //sampled_energy=(f10-f00)/ts2*rand_index+f00;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                    //1d interpolation along rand sample axis
                }
                else if(index1!=ts1-1 && index2==ts2-1)
                {
                    f00=cdf_xs_pp_tau[index1][index2];
                    f01=cdf_xs_pp_tau[index1+1][index2];
                    //sampled_energy=(f01-f00)/10.*log_energy_MeV+f00;
                    sampled_energy=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                    //1d interpolated along energy axis
                }
                else if(index1!=ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_tau[index2+1]-cdf_val_tau[index2];
                    f00=cdf_xs_pp_tau[index1][index2];
                    f01=cdf_xs_pp_tau[index1+1][index2];
                    f10=cdf_xs_pp_tau[index1][index2+1];
                    f11=cdf_xs_pp_tau[index1+1][index2+1];
                    //cout<<"f00,f01: "<<f00<<","<<f01<<endl;
                    //cout<<"f10,f11: "<<f10<<","<<f11<<endl;
                    //interpolate_sampled_energy()
                    //full 2d interpolation 
                    f_bot=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    f_top=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f10+(log_energy_MeV-log_energies_MeV[index1])/.2*f11;
                    //cout<<"f_bot,f_top: "<<f_bot<<","<<f_top<<endl;
                    //sampled_energy=(index2+1 - rand_index)/1.*f_bot+(rand_index-index2)/1.*f_top;
                    sampled_energy=(cdf_val_tau[index2+1] - rando)/d_val*f_bot+(rando-cdf_val_tau[index2])/d_val*f_top;

                    
                   
                

                   

                    //f_bot=(f01-f00)/10. * log_energy_MeV+f00;
                    //f_top=(f11-f10)/10. * log_energy_MeV+f10;
                    //sampled_energy=(f_top-f_bot)/ts2*rand_index+f_bot;
                    
                  
                    return sampled_energy;
                }

                //filler_cs=&cdf_xs_pp_tau[0][0];
            }
            case 2:
            {
                
                if(index1==ts1-1 && index2==ts2-1)
                {
                    sampled_energy=cdf_xs_pn_tau[index1][index2];
                    return sampled_energy;
                }
                else if(index1==ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_tau[index2+1]-cdf_val_tau[index2];

                    f00=cdf_xs_pn_tau[index1][index2];
                    f10=cdf_xs_pn_tau[index1][index2+1];

                    //sampled_energy=(f10-f00)/ts2*rand_index+f00;
                    sampled_energy=(cdf_val_tau[index2+1]-rando)/d_val*f00+(rando-cdf_val_tau[index2])/d_val*f10;
                    //sampled_energy=(index2+1-rand_index)/1.*f00+(rand_index-index2)/1.*f10;
                    return sampled_energy;
                    //1d interpolation along rand sample axis
                }
                else if(index1!=ts1-1 && index2==ts2-1)
                {
                    f00=cdf_xs_pn_tau[index1][index2];
                    f01=cdf_xs_pn_tau[index1+1][index2];
                    //sampled_energy=(f01-f00)/10.*log_energy_MeV+f00;
                    sampled_energy=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    return sampled_energy;
                    //1d interpolated along energy axis
                }
                else if(index1!=ts1-1 && index2!=ts2-1)
                {
                    d_val=cdf_val_tau[index2+1]-cdf_val_tau[index2];
                    f00=cdf_xs_pn_tau[index1][index2];
                    f01=cdf_xs_pn_tau[index1+1][index2];
                    f10=cdf_xs_pn_tau[index1][index2+1];
                    f11=cdf_xs_pn_tau[index1+1][index2+1];
                    //interpolate_sampled_energy()
                    //full 2d interpolation 
                    
                    //f_bot=(f01-f00)/10. * log_energy_MeV+f00;
                    //f_top=(f11-f10)/10. * log_energy_MeV+f10;
                    //sampled_energy=(f_top-f_bot)/ts2*rand_index+f_bot;

                    f_bot=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f00+(log_energy_MeV-log_energies_MeV[index1])/.2*f01;
                    f_top=(log_energies_MeV[index1+1]-log_energy_MeV)/.2*f10+(log_energy_MeV-log_energies_MeV[index1])/.2*f11;
                    sampled_energy=(cdf_val_tau[index2+1] - rando)/d_val*f_bot+(rando-cdf_val_tau[index2])/d_val*f_top;

                    //sampled_energy=(index2+1 - rand_index)/1.*f_bot+(rand_index-index2)/1.*f_top;
                    //cout<<sampled_energy<<endl;
                    return sampled_energy;
                }
                //filler_cs=&cdf_xs_pn_tau[0][0];
            }
        }
    }
    printf("failed to get sampled energy");
    exit(EXIT_FAILURE);
}

/*
int main()
{
    srand(time(NULL));
    stochastic_lepton_prop loss;
    loss.load_tables();
    double pos[3]={0,0,0};
    double traj[3]={1,0,0};
    loss.set_val(1E21,pos,traj,13,.92);
    double step=loss.get_interaction_length();
    printf("step length is %f\n",step);
    //loss.set_sto_type();
    //printf("sto type is %i",loss.sto_type);
    double e_sample=loss.get_sampled_energy();
    printf("sto type is %i and sampled energy is %f\n",loss.sto_type,e_sample);

    //cout<<loss.cdf_val_muon[10]<<endl;
    return 0;
}
*/
