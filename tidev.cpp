//**********************************************************************
//   _______     ________     
//  /_  __(_)___/ / ____/   __
//   / / / / __  / __/ | | / /
//  / / / / /_/ / /___ | |/ / 
// /_/ /_/\__,_/_____/ |___/  
// Tidal Spin Evolution
//**********************************************************************
// Copyright (C) 2013 
// Jorge Zuluaga (zuluagajorge@gmail.com, Mario Melita (melita@iafe.uba.ar)
// Pablo Cuartas (quarktas@gmail.com), Bayron Portilla (bayron@gmail.com)
//**********************************************************************
// MAIN PACKAGE FILE
//**********************************************************************
//Use: OPTIONS=-DVERBOSE make <program>.out to enable this macro
//#define VERBOSE
////////////////////////////////////////////////////////////////////////
//HEADERS
////////////////////////////////////////////////////////////////////////
//STANDARD
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <sys/time.h>

//GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_sf.h>

#include <gsl/gsl_odeiv2.h>

//LIBCONFIG
#include <libconfig.h++>

////////////////////////////////////////////////////////////////////////
//GLOBAL TYPES
////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace libconfig;
typedef void* params;
typedef double real;
typedef FILE* file;

////////////////////////////////////////////////////////////////////////
//MACROS
////////////////////////////////////////////////////////////////////////
#define setting const Setting&
#define configInit() setting CFG_ROOT=CFG.getRoot()
#define configList(var,name) setting var=CFG_ROOT[name]
#define configValue(type,var,name) type var=CFG.lookup(name)
#define configValueList(list,var,name) list.lookupValue(name,var)

//ERROR HANDLING
#define errorOff() gsl_set_error_handler_off()

//EXPONENTIATION
#define PINT(x,n) gsl_pow_int(x,n)
#define PREAL(x,y) pow(x,y)

//STRING
#define STR(string) string.c_str()

//VERBOSITY
#define Verbose printf("... ");printf

//MULTILINE STRING
#define MULTI(str) #str

////////////////////////////////////////////////////////////////////////
//CONSTANTS
////////////////////////////////////////////////////////////////////////

//NUMERICAL
//#define PI M_PIl
#define PI M_PI
#define PI2 (PI*PI)

//BEHAVIOR
#define NMAX 50

//PHYSICAL AND ASTRONOMICAL CONSTANTS
#define MSUN (GSL_CONST_MKSA_SOLAR_MASS) //kg
#define RSUN (6.96342E8) //m
#define AU (GSL_CONST_MKSA_ASTRONOMICAL_UNIT) //m
#define YEAR (365.25*GSL_CONST_MKSA_DAY) //s
#define MEARTH (5.9736E24) //kg
#define REARTH (6.371E6) //m
#define GCONST (GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT) // m^3 / (kg s^2)

#define MICRO 1E-6

//BEHAVIOR
#define MAXBODIES 10

////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////
int Option;
class Config CFG;
int NBodies,IBody;
real UL,UM,UT,GPROG;
real Gecc[11];
real TauTidal,TauTriax;

////////////////////////////////////////////////////////////////////////
//CLASSES
////////////////////////////////////////////////////////////////////////
class Body
{
public:
  string name;
  string units;

  //Physical
  real M,R;
  real rho;

  //Rheological
  real mur,alpha,tauM,tauA;
  real gapo,cosapt,sinapt;//Alpha functions 
  real A2;//A_2
  
  //Orbital parameters
  real mu;
  real a,e;
  real n,P;
  real gmt;//1.5*G*Ms^2*R^5/(a^6)

  //Momentum of Inertia: C: main, A,B: secondary
  real MoI,BmA;
  real A,B,C;
  
  //Inertia momenta
  int physicalProperties(){
    //Mean density
    if(R>0){
      rho=M/(4*PI*R*R*R);
    }else rho=0;

    //Moment of inertia
    C=MoI*M*R*R;
    A=C;
    B=C*(1.0+BmA);
    
    //Derived factor
    if(a>0){
      gmt=1.5*(mu*mu/GPROG)*R*R*R*R*R/(a*a*a*a*a*a);
    }
    else gmt=0;
  }
  
  //Orbital properties
  int orbitalProperties(){
    if(a>0){
      n=2*PI*sqrt((mu/GPROG+M)/(a*a*a));
      P=2*PI/n;
    }
    else n=P=gmt=0.0;
  }

  //Rheology properties
  int rheologyProperties(){
    gapo=gsl_sf_gamma(alpha+1);
    cosapt=cos(alpha*PI/2);
    sinapt=sin(alpha*PI/2);
    A2=(57./8)*mur/(PI*GPROG*rho*rho*R*R);
  }

  //Conversion units
  int setUnits(void);
}Bodies[MAXBODIES],Central;

////////////////////////////////////////////////////////////////////////
//ROUTINES
////////////////////////////////////////////////////////////////////////

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//I/O
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int configLoad(const char* file)
{
  if(fopen(file,"r")!=NULL) 
    CFG.readFile(file); 
  else {
    fprintf(stderr,"No configuration file '%s' found\n",file);
    exit(1);
  }
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//UTIL
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

/* ----------------------------------------------------------------------
TIME ROUTINE
---------------------------------------------------------------------- */
double Time(void)
{
  double t;
  struct timeval tiempo;
  gettimeofday(&tiempo,NULL);
  t=1e6*tiempo.tv_sec+tiempo.tv_usec;
  return t;
}

/* ----------------------------------------------------------------------
SGN FUNCTION
---------------------------------------------------------------------- */
template <typename T> 
int SGN(T val) {
  return (T(0) < val) - (val < T(0));
}

/* ----------------------------------------------------------------------
FORMAT STRING
---------------------------------------------------------------------- */
char* FRM(int n)
{
  char* cadena;
  cadena=(char*)calloc(sizeof(char),100);
  strcpy(cadena,"");
  for(int i=1;i<=n;i++){
    strcat(cadena,"%.15e ");
  }
  strcat(cadena,"\n");
  return cadena;
}

int readBodies(void)
{
  configInit();
  configList(bodies,"bodies");

  NBodies=bodies.getLength();
  IBody=1;//Object by default

  for(int i=0;i<NBodies;i++){
    //Read properties
    configValueList(bodies[i],Bodies[i].name,"name");
    configValueList(bodies[i],Bodies[i].units,"units");

    configValueList(bodies[i],Bodies[i].M,"M");
    configValueList(bodies[i],Bodies[i].R,"R");

    configValueList(bodies[i],Bodies[i].a,"a");
    configValueList(bodies[i],Bodies[i].e,"e");

    configValueList(bodies[i],Bodies[i].MoI,"MoI");
    configValueList(bodies[i],Bodies[i].BmA,"BmA");

    configValueList(bodies[i],Bodies[i].mur,"mur");
    configValueList(bodies[i],Bodies[i].alpha,"alpha");
    configValueList(bodies[i],Bodies[i].tauM,"tauM");
    configValueList(bodies[i],Bodies[i].tauA,"tauA");
    
    //Adjust units
    Bodies[i].setUnits();

    //Compute derived properties
    Bodies[i].mu=GPROG*Bodies[0].M;
    Bodies[i].physicalProperties();
    Bodies[i].orbitalProperties();
    Bodies[i].rheologyProperties();
  }

  Central=Bodies[0];

  #ifdef VERBOSE
  for(int i=0;i<NBodies;i++){
    Body b=Bodies[i];
    Verbose("Body %d: %s\n",i+1,STR(b.name));
    Verbose("\tPhysical: M = %e, R = %e\n",b.M,b.R);
    Verbose("\tConstitution: rho = %e\n",b.rho);
    Verbose("\tRheology: mur (rigidity) = %e, tau_M = %e, tau_A = %e\n",
	    b.mur,b.tauM,b.tauA);
    Verbose("\t\talpha = %e, Gamma(1+alpha) = %e\n",b.alpha,b.gapo);
    Verbose("\tInertia moment: I = (%e,%e,%e)\n",b.A,b.B,b.C);
    Verbose("\t\tgmt = %e\n",b.gmt);
    Verbose("\tOrbital basic: mu=%e, a = %e, e = %e\n",b.mu,b.a,b.e);
    Verbose("\tOrbital derived: n = %e, P = %e\n",b.n,b.P);
  }
  #endif
}

int readGEcc(string file)
{
  FILE* fh=fopen(STR(file),"r");
  real X_32[NMAX+1][NMAX+1];
  int n,m,i,j,k,q;
  real coef,tmp;
  for(j=0;j<NMAX;j++) for(k=0;k<NMAX;k++) X_32[j][k]=0.0;
  while(!feof(fh)){
    fscanf(fh,"%d %d %d %d %lf %lf",&n,&m,&j,&k,&coef,&tmp);
    if(n<-3) continue;
    else if(n==-3) if(m==2) X_32[j][k]=coef;
    else break;
  }
  for(int i=1;i<=10;i++){
    q=-3+i;
    k=2+q;
    Gecc[i]=0.0;
    for(int j=0;j<NMAX;j++){
      Gecc[i]+=X_32[j][k]*PINT(Bodies[IBody].e,j);
    }
    Gecc[i]*=Gecc[i];
  }

  #ifdef VERBOSE
  for(int i=1;i<=10;i++){
    Verbose("G[%d] = %e\n",i,Gecc[i]);
  }
  #endif

  return GSL_SUCCESS;
}

int Body::setUnits(void)
{
  real ul,um;
  if(units=="EARTH")
    {ul=(REARTH/UL);um=(MEARTH/UM);}
  else if(units=="SUN")
    {ul=(RSUN/UL);um=(MSUN/UM);}

  M*=um;
  R*=ul;
  a*=(AU/UL);
  mur=mur/(UM/(UL*UT*UT));
  tauM*=(YEAR/UT);
  tauA*=(YEAR/UT);
}

int setUnits(real ul,real um,real ut,real G)
{
  real ratio=(G/GCONST);
  if(ut==0){
    UL=ul;UM=um;GPROG=G;
    UT=sqrt(ratio*ul*ul*ul/um);
  }
  else if(ul==0){
    UM=um;UT=ut;GPROG=G;
    UL=1/pow(ratio/(um*ut*ut),1.0/3);
  }
  else if(um==0){
    UL=ul;UT=ut;GPROG=G;
    UM=ratio*ul*ul*ul/(ut*ut);
  }
  else if(G==0){
    UL=ul;UT=ut;UM=um;
    GPROG=GCONST/(ul*ul*ul/(um*ut*ut));
  }

  #ifdef VERBOSE
  Verbose("Units: UL = %e, UM = %e, UT = %e\n",UL,UM,UT);
  Verbose("Gravitational Constant: G = %e UL^3/(UM UT^2)\n",GPROG);
  #endif
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//NUMERICAL
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

/* ----------------------------------------------------------------------
SOLVE KEPLER EQUATION

Source:
Taken from http://astro.pas.rochester.edu/~aquillen/ast570/code/kepcart.cpp

supposed to be good to second order in e, from Brouwer+Clemence u0 is
first guess
---------------------------------------------------------------------- */
#define PREC_ECC_ANO 1e-14
real solveKepler(real e,real M)
{
  real du,u0,l0;
  du=1.0;
  u0=M+e*sin(M)+0.5*e*e*sin(2.0*M);
  while(fabs(du)>PREC_ECC_ANO){
    l0 = u0 - e*sin(u0);
    du = (M - l0)/(1.0 - e*cos(u0));
    u0 += du;
  }
  return u0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//PHYSICAL
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
int tidalAcceleration(double t,const double y[],double dydt[],params ps)
{
  real wterm,w220q,X220q;
  real tauTerm;
  real tauTidal,tauTriaxial,dotEtid,dota;
  Body b=Bodies[IBody];
  real n=b.n;
  int l=2,m=2,p=0,q,i;
  real chitauA;
  real R,I;

  static int ncall=1;

  #ifdef VERBOSE
  Verbose("Calling %d to fcn with t=%e, y=(%e,%e,%e)...\n",ncall++,t,y[0],y[1],y[2]);
  #endif
  
  ////////////////////////////////
  //d theta / dt = omega
  ////////////////////////////////
  dydt[0]=y[1];
  
  ////////////////////////////////
  //d omega / dt:
  ////////////////////////////////
  //Tidal torque
  tauTidal=0;
  dotEtid=0;
  wterm=(l-2*p)*n-m*y[1];
  for(int q=-2;q<=7;q++){
    i=q+3;
    w220q=wterm+q*n;
    X220q=fabs(w220q);
    chitauA=PREAL(X220q*b.tauA,b.alpha);
    R=1+b.gapo*b.cosapt/chitauA;
    I=-1.0/(X220q*b.tauM)-b.gapo*b.sinapt/chitauA;
    tauTerm=(-1.5*b.A2*I/((R+b.A2)*(R+b.A2)+I*I));
    tauTidal+=Gecc[i]*tauTerm*SGN(w220q);
    dotEtid+=Gecc[i]*tauTerm*X220q;
    #ifdef VERBOSE
    Verbose("(%d,%d): %e %e %e %e %e %e %e %e %e\n",
	   q,i,w220q,b.alpha,b.gapo,R,I,
	   tauTerm,
	   Gecc[i],
	   tauTidal,
	   dotEtid);
    #endif
  }
  tauTidal*=b.gmt;
  dotEtid*=b.gmt;
  
  //Keplerian parameters
  real M=(real)fmod((real)(n*t),(real)(2*PI));
  real E=solveKepler(b.e,M);
  real f=2*atan(sqrt((1+b.e)/(1-b.e))*sin(E/2)/cos(E/2));
  real r=b.a*(1-b.e*cos(E));

  //Triaxial torque
  tauTriaxial=1.5*(b.B-b.A)*b.mu*sin(2*(f-y[0]))/(r*r*r);

  //Store values for external use
  TauTidal=tauTidal/b.C;
  TauTriax=tauTriaxial/b.C;
  
  
  //Total acceleration
  dydt[1]=(tauTidal+tauTriaxial)/b.C;

  ////////////////////////////////
  //d E / dt:
  ////////////////////////////////
  dydt[2]=dotEtid;

  ////////////////////////////////
  //d a / dt:
  ////////////////////////////////
  dota=2*b.a*b.a/(b.mu*b.M)*dotEtid;
  dydt[3]=dota;

  ////////////////////////////////
  //d e / dt:
  ////////////////////////////////
  dydt[4]=
    1/(2*b.a)*((1-b.e*b.e)*dota/2-(1/b.M)*sqrt(b.a*(1-b.e*b.e)/b.mu)*tauTidal);

  #ifdef VERBOSE
  Verbose("KEPLER: M=%e,E=%e,f=%e,r=%e\n",M,E,f,r);
  Verbose("tauTidal=%e,tauTriax=%e,dwdt=%e,dEdt=%e\n",
	 tauTidal,tauTriaxial,dydt[1],dydt[2]);
  exit(0);
  #endif
  return GSL_SUCCESS;

}

