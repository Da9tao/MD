#include <iostream>
#include <fstream> 
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <windows.h>
using namespace std;

class MODEL
{
   public:
   MODEL();              //class constructor
   ~MODEL();            //class destructor

   //model description
   int npart;               //number of particles
   double (*x)[3];      //position of particles (in Angstroms)
   double (*v)[3];      //velocity (in A/ps)
   double (*f)[3];       //force (in kcal/mol/A)
   double (*m);          //mass (in g/mol)
   double cell[3];      //unit cell length (in A)
   bool   period;       //periodicity 0:nonperiodic 1:periodic
   double tmass;        //total mass of the system

   //MD parameters
   double dt;             //time step (in ps)
   double t;               //current time (in ps)   
   double Tset;         //Target temperature (in K)
   double Pset;         //Target pressure (in GPa)
   int nstep;              //total steps to run  
   int istep;               //current step 
   bool nvt;             //flag for NVT ensemble
   bool npt;			//flag for NPT ensemble
   
 //force field parameters
   double Do;  //in kcal/mol
   double Ro;  //in Angstrom
   double p0;  //p0= 12*Do*pow(Ro,12) in unit kcal/mol*A12
   double p1;  //p1= 12*Do*pow(Ro,6) in unit kcal/mol*A6
   double Rcut;         //cutoff distance for vdw interaction
   double rc2;            //Rcut*Rcut;
   double Ro_Rc3;
   double Ro_Rc6;
   double Ro_Rc9;
   double Ro3;

   //Nose-Hoover Thermostat
   double eta;           //damping factor in 1/ps
   double etadt;        //rate of damping in 1/ps^2
   double Q;              //mass of thermostat   in kcal/mol*ps^2
   double lns;            //size of thermostat, dimensionless 
   
   //Hoover Barostat
   double kai;          //damping factor in 1/ps
   double kaidt;        //rate of pressure damping in 1/ps^2
   double tau2;         //square relaxation time in ps^2
   
   //model properties
   double T;                 //temperature (in K)
   int df;                       //degree of freedom
   double Ek;               //kinetic energy in kcal/mol
   double Ep;               //potential energy in kcal/mol
   double Etot;             //Ek+Ep
   double Hami;              //Etot + energy of heat bath
   double P;                 //pressure    (in GPa)
   double V;                 //volume      (in A3)
   double Ptail;             //tail correction for pressure (GPa)
   double Eptail;           //tail correction for P.E. (kcal/mol)
   double strs[6];         //stress tensor in kcal/mol-A3

   //class functions
   int init();                    //initialization of variables
   int force();                //force calculation
   int integrate();          //verlet integration
   int sample();             //calculation of properties
   double myrand01(); //generate a random number between 0 and 1 [0,1]
   int dump_trj();           //output trajectory data to a file   
   int v_rescale();          //velocity rescale
   int NoseHoover();         //Nose-Hoover thermostat integrator
   ofstream outf;           //file stream of the trajectory
   ofstream logf;			//file stream of log
};


int main()
{
    MODEL model;                                 //declare a system
    LARGE_INTEGER clockFreq,startTime,endTime,delta;

    QueryPerformanceFrequency(&clockFreq); //obtain number of clock ticks per second
    cout<<"clock ticks per second: "<<clockFreq.QuadPart<<endl;
 
    model.init();                                       //initialization

    QueryPerformanceCounter(&startTime); //obtain current number of ticks count
    while(model.istep<model.nstep) {  //MD loop
        model.integrate();                         //integrate equations of motion
        if(model.istep%10==0) {
        	model.sample();                            //determine system property
		}

    }
    
	QueryPerformanceCounter(&endTime);
    delta.QuadPart = endTime.QuadPart - startTime.QuadPart;
    double ttime = (double)delta.QuadPart/clockFreq.QuadPart; //time in seconds
    printf("Total time used %.2f sec. %.2f steps per sec.",ttime,model.nstep/ttime);
    
    
    return 0;
}

int MODEL::sample()
{
    //calculation of system temperature and kinetic energy
    char null[1024];
    int i,j,k;

    //calculation of system temperature
    T = (Ek*2*4184)/(df*8.314);  //T= 2*Ek/(df*R)
      
    //Tail correction for Ep
    Eptail=(4.0/3.0)*3.1415926*npart*npart/V*Do*Ro3*(Ro_Rc9/6.0-Ro_Rc3);
    Ep+=Eptail;
    Etot = Ek + Ep;

    //Tail correction for pressure
    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
    Ptail=(8.0/3.0)*3.1415926*npart*npart/V/V*Do*Ro3*(Ro_Rc9/3.0-Ro_Rc3)*unitc;
    //calcualate pressure
    P=((Ek*2+strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
    P+= Ptail;   

    //calculate Hamiltonian
    Hami=Etot + ( df*8.314*Tset*lns/4184 + Q*eta*eta/2 )+ (Pset*V*6.0221367E2+df*8.314*Tset*kai*kai*tau2/2)/4184;

    //Display current information
    sprintf(null,"Current Step %d, Time: %.3f ps, ",istep,t);
    cout<<null;  
    sprintf(null,"T %.2f K, P %.2f MPa, V %.0f A3, ",T,P*1000,V);
    cout<<null;
    sprintf(null,"E(kcal/mol) Ek %.2f Ep %.2f Et %.4f H %.2f ",Ek,Ep,Etot,Hami);
    cout<<null;
    sprintf(null,"Etail %.0f%% Ptail %.0f%% ",Eptail/Ep*100,Ptail/P*100);
    cout<<null<<endl;
    sprintf(null,"%d %.3f %.2f %.2f %.2f %.2f %.2f %.4f %.4f",istep,t,T,P*1000,V,Ek,Ep,Etot,Hami); //pressure in MPa
    logf<<null<<endl;

    dump_trj();                     //output trajectory file
    return 0;
}

int MODEL::dump_trj()
{
    char null[1024];
    int i;    
    sprintf(null,"My MD trj: Current Step %d Time: %f ps",istep,t);
    outf<<null<<endl; //comment in trj file
    sprintf(null,"%5d",npart);
    outf<<null<<endl;    //number of particles
    for(i=0;i<npart;i++) { //position (nm) and velocity (nm/ps)
        sprintf(null,"%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f",
        1,"O","O",i+1,x[i][0]/10,x[i][1]/10,x[i][2]/10,v[i][0]/10,v[i][1]/10,v[i][2]/10);
        outf<<null<<endl;
    }
    sprintf(null,"%10.5f%10.5f%10.5f",cell[0]/10,cell[1]/10,cell[2]/10);
    outf<<null<<endl;  //box information
    return 0;
}

int MODEL::NoseHoover()
{   
    int i,k;
	
    if(npt) {
     double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
      P=((Ek*2+strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
      //Tail correction for pressure
      Ptail=(8.0/3.0)*3.1415926*npart*npart/V/V*Do*Ro3*(Ro_Rc9/3.0-Ro_Rc3)*unitc; 
      P+= Ptail;               
      kaidt = (P - Pset)*V*6.0221367E2/(npart*tau2*8.314*Tset); //in unit of 1/ps^2
      kai  += kaidt*dt; 
      V    *= exp( 3.0*kai*dt);
      cell[0]=cell[1]=cell[2]=pow(V,1.0/3.0);
    }
	
	   
    etadt =   (Ek*2-df*8.314*Tset/4184)/Q;   //in 1/ps^2
    eta    += etadt*dt;                                     //in 1/ps
    lns    += eta*dt;                                        //dimensionless
 
    for(i=0;i<npart;i++) {
        for(k=0;k<3;k++) {
           f[i][k] -= (eta+kai)*m[i]*v[i][k]/418.4;
        }
    }
    
    return 0;
}


int MODEL::force()
{  //The function determines the net force on each particle
    
    int i,j,k;
    for(i=0;i<npart;i++) {
        for(k=0;k<3;k++) f[i][k]=0;   //set forces to zero
    }
    for(i=0;i<6;i++) strs[i]=0;       //set stress to zero 

    // consider atom j at origin, the force on atom i at some position r
    // fx = - dU/dx = -(dU/dr)(dr/dx)= - (x/r)(dU/dr)
    // U = Do ( (Ro/r)^12 - 2 (Ro/r)^6) )
    // dU/dr = -12 Do/r ( (Ro/r)^12 - (Ro/r)^6) )
    // fx = 12 x Do/r^2 ( (Ro/r)^12 - (Ro/r)^6) )
    //    =  x ( 12DoRo^12/r^6 - 12DoRo^6 )/r^6 /r^2


    double r2,r2i,r6i,ff,xr[3],redu;
    Ep=0;                                             //potential energy
    for(i=0;i<npart;i++) {
        for(j=i+1;j<npart;j++) {
            for(k=0;k<3;k++) xr[k]= x[i][k] - x[j][k];  //distance vector
            if(period==1) { //periodic system, find distance within one cell  
                for(k=0;k<3;k++) { //minimum image convention
                     redu= (xr[k]/cell[k]);              //reduced coordinates
                     redu= redu - round (redu);   //between -0.5 and 0.5
                     xr[k] = redu*cell[k];               //real coordinates 
                 } 
            }                                          
            r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];   //distance squared
            if(r2<rc2 || period==0) {  //within cutoff distance            
	            r2i= 1/r2;
	            r6i= r2i*r2i*r2i;
	            ff = r6i*(p0*r6i-p1)*r2i;          //in unit of kcal/mol/A2
	            for(k=0;k<3;k++) {
	                  f[i][k]+= ff*xr[k];              //force in unit of kcal/mol/A2
	                  f[j][k]-= ff*xr[k];               //Newton's 3rd law
	            } 
                //the stress tensor 
                strs[0]+= ff*xr[0]*xr[0];  //xx in unit of kcal/mol
                strs[1]+= ff*xr[1]*xr[1];  //yy
                strs[2]+= ff*xr[2]*xr[2];  //zz
                strs[3]+= ff*xr[0]*xr[1];  //xy
                strs[4]+= ff*xr[0]*xr[2];  //xz
                strs[5]+= ff*xr[1]*xr[2];  //yz          				  
	            Ep += (p0*r6i - p1*2)*r6i/12.0;  //in unit of kcal/mol
             }
          }
    }       
    return 0;
}

int MODEL::init()
{   
    //simulation parameters     
    nstep=20000;                             //steps to run
    istep=0;                                  //current step
    npart=512;                              //number of particles
    Tset=100;                                //target temperature in K
    dt=0.01;                                 //time step in ps

    //settings for NVT ensemble
    nvt=1;                           //true for NVT simulations
    eta=etadt=lns=0;
    Q=500;                          //Nose-Hoover mass in kcal/mol*ps^2
    Hami=0;                       //Hamiltonian in kcal/mol

	//settings for NPT ensemble
    npt=0;                           //true for NPT simulations
    Pset=0.05;                    //target pressure in GPa
    kai=kaidt=0;
    tau2=1e2;                      //square relaxation time in ps^2	


    df= 3*npart;
    //allocation of memerory    
    x=new double [npart][3];        //position in Angstrom
    v=new double [npart][3];         //velocity in Angstrom/ps
    f=new double [npart][3];          //force in kcal/mol/A
    m=new double [npart];             //mass in g/mol

    //assign mass of particles
    int i,j,k;
    for(i=0;i<npart;i++) m[i]=39.948; //molecular weight of Argon        
    //One simple way to place particles in space

    //force field parameters
    Do= 0.185;    // in kcal/mol
    Ro= 3.868;    // in Angstrom
    p0= 12*Do*pow(Ro,12); //in unit kcal/mol*A12
    p1= 12*Do*pow(Ro,6);   //in unit kcal/mol*A6

    cell[0]=cell[1]=cell[2]=31.31894;      //length of unit cell in Angstroms
    period=1;                        //flag for periodicity
    Rcut=2.5*Ro;                     //cutoff distance for vdw interaction
    rc2=Rcut*Rcut;                   //Rcut2

   //calculate volume
    V=cell[0]*cell[1]*cell[2];   //in A^3
    
    //Tail correction for Ep
    Ro_Rc3=pow(Ro/Rcut,3.0);
    Ro_Rc6=Ro_Rc3*Ro_Rc3;
    Ro_Rc9=Ro_Rc6*Ro_Rc3;
    Ro3=Ro*Ro*Ro;
    Eptail=(4.0/3.0)*3.1415926*npart*npart/V*Do*Ro3*(Ro_Rc9/6.0-Ro_Rc3);

    //Tail correction for pressure
    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
    Ptail=(8.0/3.0)*3.1415926*npart*npart/V/V*Do*Ro3*(Ro_Rc9/3.0-Ro_Rc3)*unitc; 

    double sep=4;  //separation distance between two particles
    int nt=0;

    if(period) { //for periodic systems
       df -= 3; 
      //A rough method to place particles in the box
      double delta[3];
      int pps;
      pps= int(pow(npart,1.0/3.0))+1;            //particles per side
      for(k=0;k<3;k++) delta[k]=cell[k]/pps; //spacing of particles
      for(i=0;i<pps;i++) {
          for(j=0;j<pps;j++) {
              for(k=0;k<pps;k++) {
                  x[nt][0]=i*delta[0];
                  x[nt][1]=j*delta[1];
                  x[nt][2]=k*delta[2];
                  nt++;                                       //number of particles placed
                  if(nt==npart) i=j=k=pps;        //nt has reached npart, exit loops
              }
          }
      }        	
	} else {  //nonperiodic systems
	    for(i=0;i<npart;i++) {
	        for(j=0;j<=i;j++) {
	            for(k=0;k<=j;k++) {
	                x[nt][0]=i*sep;
	                x[nt][1]=j*sep;
	                x[nt][2]=k*sep;
	                nt++;
	                if(nt==npart) i=j=k=npart;  //stop when nt==npart
	            }
	        }
	    }
	}    
	//Assign velocities
    double cmv[3],sumv2,Ti,fs;
    srand(123);                                    //seed for random number
    cmv[0]=cmv[1]=cmv[2]=sumv2=tmass=0;
                                        //degree of freedom

    for(i=0;i<npart;i++) {
        tmass += m[i];                           //total mass
        for(k=0;k<3;k++) {
            v[i][k]= myrand01()-0.5;        //random number between -0.5 and 0.5
            cmv[k]+= m[i]*v[i][k];            //center of mass velocity
            sumv2 += m[i]*v[i][k]*v[i][k];// m*v2
        }
    }
    
    for(k=0;k<3;k++) { cmv[k]/=tmass; sumv2 -= tmass*cmv[k]*cmv[k]; }
    Ti = sumv2/(df*8.314*0.1);            //initial temperature from random vel
    fs = sqrt(Tset/Ti);                           //scale factor
    
    sumv2 =0; 
    for(i=0;i<npart;i++) {
        for(k=0;k<3;k++) {
            v[i][k]  = (v[i][k]-cmv[k]) *fs;  //rescale initial vel, remove C.M. vel.
            sumv2 += m[i]*v[i][k]*v[i][k];// m*v2
        }
    }
    Ek=sumv2/(2*418.4);      //in kcal/mol;
    //calculation of system temperature
    T = (Ek*2*4184)/(df*8.314);  //T= 2*Ek/(df*R)

    outf.open("mytrj.gro",ios::out);   //output filename for trajectory
    logf.open("mylog.txt",ios::out);   //log file
    force();                                            //initial force calculation 
    
    return 0;
}

int MODEL::v_rescale()
{
	int nn=0.1*npart;
    int k,i;
	double cmv[3],tmass;
	double sumv2,Ti,fs,df;
    int id[nn];
 
 	tmass=sumv2=cmv[0]=cmv[1]=cmv[2]=0;
	df=3*nn;
    for(i=0;i<nn;i++) {                           
	    tmass += m[i]; //total mass
        for(k=0;k<3;k++) {
            cmv[k]+= m[i]*v[i][k];            //center of mass velocity
            sumv2 += m[i]*v[i][k]*v[i][k];// m*v2
        }
    }
	 
    for(k=0;k<3;k++) { cmv[k]/=tmass; sumv2 -= tmass*cmv[k]*cmv[k]; }
    Ti = sumv2/(df*8.314*0.1);            //current temperature from vel
    fs = sqrt(Tset/Ti);                           //scale factor


    for(i=0;i<nn;i++) {
        for(k=0;k<3;k++) {
            v[i][k]  = (v[i][k]-cmv[k]) *fs;  //rescale initial vel, remove C.M. vel.
            v[i][k] += cmv[k];
        }
    }

	return 0;
} 

int MODEL::integrate()
{    //velocity Verlet integration for particle position and velocity   
    int i,k;
    double tmp;
    
    for(i=0;i<npart;i++) {
    	tmp = dt*418.4/(2*m[i]);
        for(k=0;k<3;k++) {
           // force in kcal/mol/A, mass in g/mol, x in A, dt in ps, v in A/ps
           x[i][k] +=  (v[i][k] + f[i][k]*tmp)*dt;
           v[i][k] +=  f[i][k]*tmp;   //contribution to vel from f[i]
        }
    }
    force();                                         //cal forces at r(t+dt)
    if(nvt||npt) NoseHoover();
    
    Ek=0;
    
    for(i=0;i<npart;i++) {
    	tmp = dt*418.4/(2*m[i]);
        for(k=0;k<3;k++) {
           v[i][k] += f[i][k]*tmp; //contribution to vel from fm[i]
           Ek += m[i]*v[i][k]*v[i][k];
        }
    }
    Ek /= (2*418.4);        //in kcal/mol
 
	//v_rescale(); //rescale velocity to set temperature        
    istep++;                                         //current step        
    t=istep*dt;                                     //current time in ps
    return 0;
}

MODEL::MODEL()
{
    
};

MODEL::~MODEL()
{
    delete [] x;
    delete [] v;    
    delete [] f;    
    delete [] m;
    outf.close();
    logf.close();
};

double MODEL::myrand01()
{
    return rand()/double(RAND_MAX);  
     //returns a number between 0 (inclusive) and 1 (inclusive)
}
