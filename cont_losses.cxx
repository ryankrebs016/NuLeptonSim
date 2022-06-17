//class for continuous lepton losses during propagation
#include<stdio.h>
#include<iostream>
#include<math.h>
#include"Constantes.hh"


using namespace std; 

class continuous_loss_prop 
{
    public:
    double energy_GeV;
    double dens;
    int p_type;
    const double m_frac = 1e-2;
    const double t_frac=1e-2;
    int loss_model;
    double dL;
    double elost;

    double calc_e_loss();
    double get_interaction_length();
    double get_energy_loss();
    double funcalph(int lyr);
    double beta9fit(int lyr);
    void set_values(double energy,double t_dens,int type,int t_loss_model);
    double delta(double en);
    void set_zero();
};

void continuous_loss_prop::set_zero()
{
    dL=0;
    elost=0;
}
double continuous_loss_prop::get_interaction_length()
{
    calc_e_loss();

    if(p_type==13)
    {
        dL=energy_GeV/(dens*elost)*m_frac;
    }
    else if(p_type==15)
    {
        dL=energy_GeV/(dens*elost)*t_frac;
        
    }
    return dL;

}
double continuous_loss_prop::get_energy_loss()
{
    
    
    return dL*dens*elost;
    

}

double continuous_loss_prop::calc_e_loss()
{
  double f;
  //double z = 0.;
  //  dE/dX =      E*beta(E)      +     alpha(E)
  
  // this is a super-kludgy way to account for beta being different for iron <A>=56.84, <Z>=26; rock <A>=22, <Z>=11; and water <A>=11.9, <Z>=6.6
  // material properties are not tracked in this simulation but densities are.

  int lyr = 0; // initialize to iron
  if (dens < 7.75) lyr = 1; // density jump between outer core and mantle
  if (dens < 2.0) lyr = 2;  // density jump between rock and water 
  if(p_type==13)
  {
  double factor[3] = {0.9304, 1.0, 1.1092}; // ratio Z/A for iron, rock, and water divided by Z/A=0.5 for rock   

  // The correction below was based on the claim that the photonuclear energy loss is proportional to <A> as well as the density in Palomares-Ruiz, Irimia, & Weiler, Phys. Rev. D, 73, 083003 (2006)
  // Searching through the references, this claim is demonstrably false. See S. I. Dutta, M. H. Reno, I. Sarcevic, and D. Seckel, Phys. Rev. D 63, 094020 (2001).
  // Earlier runs of the code used the line below but it has been commented out.
  // if(dens<1.1) factor = 0.55; // This kluge corrects for the change of <A>=22 in rocks vs <A>=12 of H2O
  //printf("factor %1.2f \n", factor);
  
  //f = E * beta9fit(&E,&lyr,ELOSSmode) + mfuncalph(&E, &lyr,type);
  f = 2e-3*factor[lyr] + energy_GeV * beta9fit(lyr);
  
  
	//f = E*(emlost->Eval(E)) + funca->Eval(E);
  // cout << "\tE " << E << endl;
  //cout << "\tlyr " << lyr << " dens = " << dens << endl;
  // cout << " f = " << f << endl;
  // cout << " E * beta9fit(&E,&lyr) + mfuncalph(&E, &lyr); " << E * beta9fit(&E,&lyr) + mfuncalph(&E, &lyr) << endl << endl;
  // cout << " mfuncalph(&E,&lyr) " << E << " " << lyr << " " << mfuncalph(&E, &zlyr << endl << endl;
  //cout << " beta9fit(&E,&lyr) " << E << " " << beta9fit(&E,&lyr,0) << " " << beta9fit(&E,&lyr,1) << endl << endl;
  }
  if(p_type==15)
  {
    f = energy_GeV * beta9fit(lyr) + funcalph(lyr);

    //f = E*(emlost->Eval(E)) + funca->Eval(E);
    // cout << "\tE " << E << endl;
    //cout << "\tlyr " << lyr << " dens = " << dens << endl;
    // cout << " f = " << f << endl;
    // cout << " E * beta9fit(&E,&lyr) + funcalph(&E, &lyr); " << E * beta9fit(&E,&lyr) + funcalph(&E, &lyr) << endl << endl;
    // cout << " funcalph(&E,&lyr) " << E << " " << lyr << " " << funcalph(&E, &zlyr << endl << endl;
    //cout << " beta9fit(&E,&lyr) " << E << " " << beta9fit(&E,&lyr,0) << " " << beta9fit(&E,&lyr,1) << endl << endl;
  }
  //if(f==0){cout<<"elost is zero"<<endl;}
  elost=f;
  return f;

}
void continuous_loss_prop::set_values(double energy,double t_dens,int type,int t_loss_model)
{
    energy_GeV=energy;
    dens=t_dens;
    p_type=type;
    loss_model=t_loss_model;
    
}

double continuous_loss_prop::funcalph(int lyr)
{
  double f;
  double p,b,b2,gamma,EE,X;
  double factor[3]={0.9304,1.0,1.1092};

  if(p_type==15)
  {
      p=sqrt(energy_GeV*energy_GeV-mtau2);
      b=p/energy_GeV;
      b2=b*b;
      gamma=energy_GeV/mtau;
      EE=Cbb2*p*p/(me2+mtau2+Cbb2*energy_GeV);
      X=log10(b*gamma);
      
  }
  else if( p_type==13)
  {
    
  
      p=sqrt(energy_GeV*energy_GeV-mmuon2);
      b=p/energy_GeV;
      b2=b*b;
      gamma=energy_GeV/mmuon;
      EE=Cbb2*p*p/(me2+mmuon2+Cbb2*energy_GeV);
      X=log10(b*gamma);
      
  
  }
  //printf("p,b,b2,gamma,EE,X %f,%f,%f,%f,%f,%f\n",p,b,b2,gamma,EE,X);

  f=Cbb1/(b2)*(log(Cbb2*b2*gamma*gamma/I2)-2.*b2+EE*EE/(4.*energy_GeV*energy_GeV)-delta(X));
  
  f*=factor[lyr];
  return f;

}
double continuous_loss_prop::delta(double en)
{
  double f=0;
  
  if(en > X1)
  {
    f=4.6052*en+CC;
  }
  else if((en>X0)&&(en<X1))
  {
    f=4.6052*en+CC+aa*pow((X1-en),mm);
  }
  
  return f;
}

double continuous_loss_prop::beta9fit(int lyr)
{
    double f=0.;
    double b0 = 0.;
    double b1 = 0.;
    double b2 = 0.;
    double b3 = 0.;

    if(p_type==15||p_type==16) //tau
    {
      
        
        /* ALLM */
        if(loss_model==0)
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
        if(loss_model==1)
        {
        b0=-4.77043758142e-08;
        b1=1.9031520827e-07;
        b2=0.0469916563971;
        b0 = -8.61998129e+00;
        b1 =  1.57820040e-01; 
        b2 = -2.24340096e-03;
        //printf("ASW \n");
        }
        double log10E = log10(energy_GeV);
        double f_phot = pow(10., b0 + b1*(log10E+9.) + b2*(log10E+9.)*(log10E+9.));
        //f=b0+b1*pow(x[0],b2);
        //printf("%1.2e \n", f);
        
        double f_brem = pBrem[lyr][0]*pow(1.-exp(-pow((log10E)/pBrem[lyr][1], pBrem[lyr][2])),pBrem[lyr][3]);
        double f_pair = pPair[lyr][0]*pow(1.-exp(-pow((log10E)/pPair[lyr][1], pPair[lyr][2])),pPair[lyr][3]);

        // cout << "\tpar[0] " << par[0] << endl;
        // cout << "\tBrem " << pBrem[par[0]][0] << " " <<  pBrem[par[0]][1] << " " <<  pBrem[par[0]][2] << " " <<  pBrem[par[0]][3] << endl;
        // cout << "\tPair " << pPair[par[0]][0] << " " <<  pPair[par[0]][1] << " " <<  pPair[par[0]][2] << " " <<  pPair[par[0]][3] << endl;
        // cout << "\ttest " << 1.-exp(pow(log10(x[0])/pBrem[par[0]][1], pBrem[par[0]][2])) << endl;
        // cout << "\tpair, brem, photoN " << f_pair << " " << f_brem << " " << f_phot << endl;
        f=f_pair + f_brem + f_phot;
        

    };
    if(p_type==13||p_type==14) //muon 
    {

            /*Energy loss parameters here are the sum ofbremmstrahlung, pair production, and photonuclear interactions*/
            /*Photonuclear losses are characterized using the 2 models below*/
           

            //int lyr = par[0];
                
            /* BB */
            /* Bezrukov and Bugaev Model for photonuclear losses*/
            /* L. B. Bezrukov and E. V. Bugaev, Sov. J. Nucl. Phys. 33, 635 (1981). */
            if(loss_model==1)
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
            if(loss_model==0)
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
            double log10E = log10(energy_GeV);
            f = b0+b1*log10E+b2*log10E*log10E+b3*log10E*log10E*log10E;

            
    }
    
    return f;
}

/*
how to use it in the main sim.
int main()
{
    continuous_loss_prop cont;
    int loop=100000;
    double energy=1E19;
    cont.set_values(1E19,.92,15,0);
   
    int loop_num=0;
    double dist=0;
    for(int i=0;i<loop;i++)
    {
        loop_num++;
        dist+=cont.get_interaction_length();
        
        energy=energy-cont.get_energy_loss();
        //cout<<cont.dL<<","<<cont.get_energy_loss()<<endl;;
        cout<<dist<<","<<energy<<endl;
        cont.set_values(energy,.92,15,0);
        if(energy<1E16) break;

    }
    cout<<"traveled "
    cout<<loop_num<<"loops"<<endl;
   return 0;

}
*/
