// **********************************************************************
//   _______     ________     
//  /_  __(_)___/ / ____/   __
//   / / / / __  / __/ | | / /
//  / / / / /_/ / /___ | |/ / 
// /_/ /_/\__,_/_____/ |___/  
// Tidal Spin Evolution
// **********************************************************************
// Copyright (C) 2013 
// Jorge Zuluaga (zuluagajorge@gmail.com, Mario Melita (melita@iafe.uba.ar)
// Pablo Cuartas (quarktas@gmail.com), Bayron Portilla (bayron@gmail.com)
// **********************************************************************
// MAIN PACKAGE FILE
// **********************************************************************
// Use: OPTIONS=-DVERBOSE make <program>.out to enable this macro
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
#include <gsl/gsl_integration.h>
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

//ALPHA
#define ALPHA(a1,a2) (a1<=a2?a1/a2:a2/a1)

////////////////////////////////////////////////////////////////////////
//CONSTANTS
////////////////////////////////////////////////////////////////////////

//NUMBER OF VARIABLES
#define NUMVARS 8 //5 elements + theta + omega + Etid

//METHOD OF INTEGRATION
#define ODEMETHOD gsl_odeiv2_step_rk4 

//NUMERICAL
//#define PI M_PIl
#define PI M_PI
#define PI2 (PI*PI)
#define D2R PI/180

//BEHAVIOR
#define NMAX 50

//PHYSICAL AND ASTRONOMICAL CONSTANTS
#define MSUN (GSL_CONST_MKSA_SOLAR_MASS) //kg
#define RSUN (6.96342E8) //m
#define AU (GSL_CONST_MKSA_ASTRONOMICAL_UNIT) //m
#define HOURS (3600.0) //s
#define YEAR (365.25*GSL_CONST_MKSA_DAY) //s
#define MEARTH (5.9736E24) //kg
#define REARTH (6.371E6) //m
#define GCONST (GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT) // m^3 / (kg s^2)

#define MICRO 1E-6

//BEHAVIOR
#define MAXBODIES 10

//MODES
#define MODE_TIDAL 0
#define MODE_SECULAR 1
#define MODE_TIDALSECULAR 2

////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////
int Mode;
int Option;
class Config CFG;
int NBodies,NPlanets,IBody;
real UL,UM,UT,GPROG;
real Gecc[11];
real TauTidal,TauTriax;

////////////////////////////////////////////////////////////////////////
//CLASSES
////////////////////////////////////////////////////////////////////////
class Body
{
public:
  int active;
  int tidal;

  string name;
  string units;

  //Physical
  real M,R;
  real rho;

  //Rheological
  real mur,alpha,tauM,tauA;
  real thetaini,Pini;
  real gapo,cosapt,sinapt;//Alpha functions 
  real A2;//A_2
  
  //Orbital parameters
  real mu;
  real a,e;
  real n,P;
  real I,Om,w;
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
PRINT VECTOR
---------------------------------------------------------------------- */
void fprintf_vec(FILE* stream,const char* fmt,const double x[],int end,int ini=0)
{
  for(int i=ini;i<end;i++){
    fprintf(stream,"%d:",i);
    fprintf(stream,fmt,x[i]);
  }
  fprintf(stream,"\n");
}

/* ----------------------------------------------------------------------
TIME ROUTINE
---------------------------------------------------------------------- */
double** matrixAlloc(int n,int m)
{
  int i;
  double **M;
  M=(double**)calloc(n,sizeof(double*));
  for(i=0;i<n;i++)
    M[i]=(double*)calloc(m,sizeof(double));
  return M;
}

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
  NPlanets=NBodies-1;
  IBody=1;//Object by default

  for(int i=0;i<NBodies;i++){
    //Read properties
    configValueList(bodies[i],Bodies[i].active,"active");
    configValueList(bodies[i],Bodies[i].tidal,"tidal");

    configValueList(bodies[i],Bodies[i].name,"name");
    configValueList(bodies[i],Bodies[i].units,"units");

    configValueList(bodies[i],Bodies[i].M,"M");
    configValueList(bodies[i],Bodies[i].R,"R");

    configValueList(bodies[i],Bodies[i].a,"a");
    configValueList(bodies[i],Bodies[i].e,"e");
    configValueList(bodies[i],Bodies[i].I,"I");
    configValueList(bodies[i],Bodies[i].Om,"Om");
    configValueList(bodies[i],Bodies[i].w,"w");

    configValueList(bodies[i],Bodies[i].MoI,"MoI");
    configValueList(bodies[i],Bodies[i].BmA,"BmA");

    configValueList(bodies[i],Bodies[i].mur,"mur");
    configValueList(bodies[i],Bodies[i].alpha,"alpha");
    configValueList(bodies[i],Bodies[i].tauM,"tauM");
    configValueList(bodies[i],Bodies[i].tauA,"tauA");

    configValueList(bodies[i],Bodies[i].thetaini,"thetaini");
    configValueList(bodies[i],Bodies[i].Pini,"Pini");
    
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
    Verbose("\tOrbital basic: mu=%e, a = %e, e = %e, i = %f deg, xi = %f deg, Om = %f deg, w = %f deg\n",
	    b.mu,b.a,b.e,b.i,b.xi,b.Om,b.w);
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
  thetaini*=(D2R);
  Pini*=(HOURS/UT);
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
////////////////////////////////////////////////////////////////////////
//LAPLACE COEFFICIENTS
////////////////////////////////////////////////////////////////////////
#define JEFF(j) ((int)(j+4))
#define SEFF(s) ((int)(s-0.5))
double coefLaplace(double th,void *param)
{
  double* ps=(double*) param;
  double j=ps[0];
  double s=ps[1];
  double alpha=ps[2];

  double f=(1/pow((1.0-2.0*alpha*cos(th)+alpha*alpha),s))*cos(j*th);
  
  return f;
}
class LaplaceCoefficients
{
public:
  double Alpha;
  double Mp;
  double time;
  gsl_matrix **Bjs;
  gsl_integration_workspace* IntWork;
  gsl_function Func;
  double Laplparam[3];
  
  double df;
  double df0;
  double df1, df2, df3, df4;
  double df01, df02;
  double df20, df22;
  double df10, df11, df12;
  double df13, df14;
  double df31, df32, df31b, df32b;
  double df100, df101, df102, df103, df104;

  int set(double alpha,double mp)
  {
    Mp=mp;
    Alpha=alpha;
    #ifdef VERBOSE
    Verbose("Initial alpha:%e\n",Alpha);
    #endif
    Bjs=(gsl_matrix**)calloc(5,sizeof(gsl_matrix*));
    IntWork=gsl_integration_workspace_alloc(1000);
    Func.function=&coefLaplace;
    Func.params=&Laplparam;
    update(alpha);
    return 0;
  }

  int update(double alpha)
  {
    Alpha=alpha;
    for(int i=0;i<=4;i++){
      Bjs[i]=gsl_matrix_alloc(11,11);
      gsl_matrix_set_all(Bjs[i],2E100);
    }
  }

  double Bsget(int order,int j,double s)
  {
    return gsl_matrix_get(Bjs[order],JEFF(j),SEFF(s));
  }

  double Bsset(int order,int j,double s,double f)
  {
    gsl_matrix_set(Bjs[order],JEFF(j),SEFF(s),f);
  }


  double bsget(int order,int j,double s,double alpha)
  {
    double f;
    f=Bsget(order,j,s);
    f>1E100&&bsfunc(order,j,s,alpha,&f);
    return f;
  }

  int bsfunc(int order,int j,double s,double alpha,double *f)        
  {
    switch(order){
    case 0:{
      double integral,error;
      Laplparam[0]=j;
      Laplparam[1]=s;
      Laplparam[2]=alpha;
      gsl_integration_qags(&Func,0.0,PI,0.0,1E-6,1000,IntWork,&integral,&error);
      *f=(2/PI)*integral;
      break;
    }
    case 1:{
      double f1,f2,f3;
      f1=bsget(0,j,s+1.0,alpha);
      f2=bsget(0,j-1.0,s+1.0,alpha);
      f3=bsget(0,j+1.0,s+1.0,alpha);
      *f=s*(f2-2.0*alpha*f1+f3);
      break;
    }
    case 2:{
      double f0=bsget(0,j,s+1.0,alpha);
      double f1=bsget(1,j,s+1.0,alpha);
      double f2=bsget(1,j-1.0,s+1.0,alpha);
      double f3=bsget(1,j+1.0,s+1.0,alpha);
      *f=s*(f2-2.0*alpha*f1-2.0*f0+f3);
      break;
    }
    case 3:{
      double f0=bsget(1,j,s+1.0,alpha);
      double f1=bsget(2,j,s+1.0,alpha);
      double f2=bsget(2,j-1.0,s+1.0,alpha);
      double f3=bsget(2,j+1.0,s+1.0,alpha);
      *f=s*(f2-4.0*f0-2.0*alpha*f1+f3);
      break;
    }
    case 4:{
      double f0=bsget(2,j,s+1.0,alpha);
      double f1=bsget(3,j,s+1.0,alpha);
      double f2=bsget(3,j-1.0,s+1.0,alpha);
      double f3=bsget(3,j+1.0,s+1.0,alpha);
      *f=s*(f2-6.0*f0-2.0*alpha*f1+f3);
      break;
    }
    }
    Bsset(order,j,s,*f);
  }

  void coefcs()
  {
    double t1,t2;
    double alpha=Alpha;

    t1=Time();
    df=bsget(0,0.0,0.5,alpha);
    df0=bsget(0,1.0,1.5,alpha);
    df01=bsget(1,1.0,1.5,alpha);
    df02=bsget(2,1.0,1.5,alpha);

    df1=bsget(1,0.0,0.5,alpha);
    df2=bsget(2,0.0,0.5,alpha);
    df3=bsget(3,0.0,0.5,alpha);
    df4=bsget(4,0.0,0.5,alpha);      

    df20=bsget(0,0.0,2.5,alpha);
    df22=bsget(0,2.0,2.5,alpha);

    df10=bsget(0,1.0,0.5,alpha);
    df11=bsget(1,1.0,0.5,alpha); 
    df12=bsget(2,1.0,0.5,alpha);
    df13=bsget(3,1.0,0.5,alpha);
    df14=bsget(4,1.0,0.5,alpha);
    
    df31=bsget(1,0.0,1.5,alpha);
    df32=bsget(2,0.0,1.5,alpha);
    df31b=bsget(1,2.0,1.5,alpha);
    df32b=bsget(2,2.0,1.5,alpha);
    
    df100=bsget(0,2.0,0.5,alpha);
    df101=bsget(1,2.0,0.5,alpha);
    df102=bsget(2,2.0,0.5,alpha);
    df103=bsget(3,2.0,0.5,alpha);
    df104=bsget(4,2.0,0.5,alpha);

    t2=Time();
    time=(t2-t1)*1E-6;

    #ifdef VERBOSE
    Verbose("Alpha = %e\n",alpha);
    Verbose("df,df0,df01,df02 = %e %e %e %e\n",df,df0,df01,df02);
    Verbose("df1,df2,df3,df4 = %e %e %e %e\n",df1,df2,df3,df4);
    Verbose("df20,df22 = %e %e\n",df20,df22);
    Verbose("df10,df11,df12,df13,df14 = %e %e %e %e %e\n",df10,df11,df12,df13,df14);
    Verbose("df31,df32,df31b,df32b = %e %e %e %e\n",df31,df32,df31b,df32b);
    Verbose("df100,df101,df102,df103,df104 = %e %e %e %e %e\n",df100,df101,df102,df103,df104);
    #endif

  }
};
#define LC LaplaceCoefficients
LC** laplaceAlloc(int n,int m)
{
  int i;
  LC** M;
  M=(LC**)calloc(n,sizeof(LC*));
  for(i=0;i<n;i++)
    M[i]=(LC*)calloc(m,sizeof(LC));
  return M;
}

////////////////////////////////////////////////////////////////////////
//SECULAR VARIATION OF ELEMENTS
////////////////////////////////////////////////////////////////////////
class SecularEvolution
{
public:
  int Np;
  int Nptid;
  double* X;
  int* Iplanets;
  int* Iplanetstid;
  LaplaceCoefficients **Laplaces;
  LaplaceCoefficients *laplace;

  int set(int np,int iplanets[],int nptid,int iplanetstid[],
	  double *x,LaplaceCoefficients** ls)
  {
    Np=np;
    Iplanets=iplanets;
    
    Nptid=nptid;
    Iplanetstid=iplanetstid;
    
    X=(double*)calloc(NUMVARS*Np,sizeof(double));
    for(int i=0;i<NUMVARS*Np;i++) X[i]=x[i];
    Laplaces=ls;

    #ifdef VERBOSE
    Verbose("\tSecular evolution properties:\n");
    Verbose("\tNumber of planets: %d\n",Np);
    Verbose("\tInitial elements: ");
    fprintf_vec(stdout,"%e ",X,NUMVARS*Np);
    #endif
    return 0;
  }

  int update(const double *x)
  {
    for(int i=0;i<NUMVARS*Np;i++){X[i]=x[i];}
  }
  
  double f2(double alpha){
    return alpha*laplace->df0/8.0;
  }

  double f3(double alpha){
    return -0.5*alpha*laplace->df0;
  }

  double f4(double alpha){
    return (4.0*(alpha*alpha*alpha)*laplace->df3 + (alpha*alpha*alpha*alpha)*laplace->df4)/128.0;
  }

  double f5(double alpha){
    return (4.0*alpha*laplace->df1 + 14.0*(alpha*alpha)*laplace->df2 +
	    8.0*(alpha*alpha*alpha)*laplace->df3 + (alpha*alpha*alpha*alpha)*laplace->df4)/32.0;
  }

  double f6(double alpha){
    return (24.0*alpha*laplace->df1 + 36.0*(alpha*alpha)*laplace->df2 +
	    12.0*(alpha*alpha*alpha)*laplace->df3 + (alpha*alpha*alpha*alpha)*laplace->df4)/128.0;
  }

  double f7(double alpha){
    return (-2.0*alpha*laplace->df0 - 4.0*(alpha*alpha)*laplace->df01 -
	    (alpha*alpha*alpha)*laplace->df02)/8.0;
  }

  double f8(double alpha){
    return (3.0/8.0)*(alpha*alpha)*laplace->df22 +
      (3.0/4.0)*(alpha*alpha)*laplace->df20;
  }

  double f9(double alpha){
    return 0.5*alpha*laplace->df0 + (3.0/4.0)*(alpha*alpha)*laplace->df22 +
      (15.0/4.0)*(alpha*alpha)*laplace->df20;
  }

  double f10(double alpha){
    return 0.25*(2.0*laplace->df10 - 2.0*alpha*laplace->df11 -
		 (alpha*alpha)*laplace->df12);
  }

  double f11(double alpha){
    return (-4.0*(alpha*alpha)*laplace->df12 - 6.0*(alpha*alpha*alpha)*laplace->df13 -
	    (alpha*alpha*alpha*alpha)*laplace->df14)/32.0;
  }

  double f12(double alpha){
    return (4.0*laplace->df10 - 4.0*alpha*laplace->df11 - 22.0*(alpha*alpha)*laplace->df12 -
	    10.0*(alpha*alpha*alpha)*laplace->df13 - (alpha*alpha*alpha*alpha)*laplace->df14)/32.0;
  }

  double f13(double alpha){
    return 0.5*(alpha*alpha)*(laplace->df31 + laplace->df31b) +
      (alpha*alpha*alpha)*(laplace->df32 + laplace->df32b)/8.0;
  }

  double f14(double alpha){
    return alpha*laplace->df0;
  }

  double f15(double alpha){
    return 0.5*alpha*laplace->df0 + (alpha*alpha)*laplace->df01 +
      0.25*(alpha*alpha*alpha)*laplace->df02;
  }

  double f16(double alpha){
    return -0.5*alpha*laplace->df0 - 3.0*(alpha*alpha)*laplace->df20 -
      1.5*(alpha*alpha)*laplace->df22;
  }

  double f17(double alpha){
    return (12.0*laplace->df100 - 12.0*alpha*laplace->df101 + 
	    6.0*(alpha*alpha)*laplace->df102 +
	    8.0*(alpha*alpha*alpha)*laplace->df103 + (alpha*alpha*alpha*alpha)*laplace->df104)/64.0;
  }

  double f18(double alpha){
    return 0.75*alpha*laplace->df0 + 0.5*(alpha*alpha)*laplace->df01 +
      (alpha*alpha*alpha)*laplace->df02/16.0;
  }

  double f19(double alpha){
    return -0.5*(alpha*alpha)*laplace->df31 - (alpha*alpha*alpha)*laplace->df32/8.0;
  }

  double f20(double alpha){
    return (alpha*alpha*alpha)*laplace->df02/16.0;
  }

  double f21(double alpha){
    return -1.5*alpha*laplace->df0 - (alpha*alpha)*laplace->df01 -
      (alpha*alpha*alpha)*laplace->df02/8.0;
  }

  double f22(double alpha){
    return -(alpha*alpha)*laplace->df31 - (alpha*alpha*alpha)*laplace->df32/4.0;
  }

  double f23(double alpha){
    return -(alpha*alpha )*laplace->df31b - (alpha*alpha*alpha)*laplace->df32b/4.0;
  }

  double f24(double alpha){
    return (alpha*alpha)*laplace->df31 + (alpha*alpha*alpha)*laplace->df32/4.0;
  }

  double f25(double alpha){
    return - (alpha*alpha*alpha)*laplace->df02/8.0;
  }

  double f26(double alpha){
    return 0.5*alpha*laplace->df0 + 0.75*(alpha*alpha)*laplace->df20 +
      1.5*(alpha*alpha)*laplace->df22;
  }

  int secular(double *dxdt)
  {
    /*
    //Pert. INTERNA O EXTERNA 
    //d/dt de e, i, \varpi, \Omega 
    //4o orden en e / si = sin((1/2)*xi)           
    */
    double a,e,xi,Om,w;
    double aM,eM,xiM,OmM,wM,mp;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //VARIABLE DECLARATION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    int i,ki,kk,ikk,j,kj;
    double aux,si,sini,Omega;
    double dxdtp1,dxdtp2,dxdtp3,dxdtp4;
    double wp0;
    double Op0;
    int nkk;
    double wp[16],Op[16],xep[16],xip[16],Rcos[16],Rsin[16];
    double mpla;
    double apla;
    double epla;
    double sipla;

    double alpha;

    double Cn;
    double Cwe;
    double Cwi;
    double Ce;
    double Ciw;
    double Cin;

    double nodo_pla;
    double w_pla;
    double nodo;
    double dc0;
    double c0;
    double c1;
    double dc1;
    double a1;
    double c2;
    double dc2;
    double a2;
    double c3;
    double a3;
    double c4;
    double dc4;
    double a4;
    double c5;
    double dc5;
    double a5;
    double dc6;
    double a6;
    double c7;
    double dc7;
    double a7;
    double c8;
    double dc8;
    double a8;
    double c9;
    double dc9;
    double a9;
    double c10;
    double dc10;
    double a10;
    double dc11;
    double a11;
    double c12;
    double a12;
    double c13;
    double a13;
    double dc15;
    double a15;
    double facint=1;

    for(i=0;i<Np;i++){
      ki=NUMVARS*i;
      //fprintf(stdout,"Perturbed Planet %d (ki=%d):\n",i,ki);
      a=X[0+ki];
      e=X[1+ki];
      xi=X[2+ki];
      Om=X[3+ki];
      w=X[4+ki];

      //Set to zero xi,Omega,w+Omega
      /*
      dxdt[2+ki]=0.0;
      dxdt[3+ki]=0.0;
      dxdt[4+ki]=0.0;
      */

      aux = sqrt(1.0-e*e);
      si = sin(xi/2.0);
      sini = sin(xi);
      Omega = sqrt(GPROG*Bodies[0].M/(a*a*a));

      for(j=0;j<Np;j++){
	if(i==j) continue;
	kj=NUMVARS*j;
	//fprintf(stdout,"\tPerturber Planet %d (ki=%d):\n",j,kj);
	aM=X[0+kj];
	eM=X[1+kj];
	xiM=X[2+kj];
	OmM=X[3+kj];
	wM=X[4+kj];
	mp=Laplaces[i][j].Mp;
	//fprintf(stdout,"\tmp =%e\n",mp);

	laplace=&Laplaces[i][j];

	wp0 = 0.0;
	Op0 = 0.0;
	nkk = 15;

	for(kk=1;kk<=nkk;kk++){
	  wp[kk] = 0.0;
	  Op[kk]= 0.0;
	  xep[kk]= 0.0;
	  xip[kk] = 0.0;
	  Rcos[kk] = 0.0;
	  Rsin[kk] = 0.0;
	}

	mpla = GPROG*mp;
	apla = aM;
	epla = eM;
	sipla = sin(xiM/2.0);

	alpha=ALPHA(a,apla);
	/*
	alpha=a/apla;
	facint=1.0;
	if(a>apla){
	  alpha=1.0/alpha;
	  //facint*=alpha;
	  //mpla*=alpha;
	  #ifdef VERBOSE
	  fprintf(stdout,"Pertubed %d (a=%e), Perturber %d (a=%e), alpha = %e, facint = %e\n",
		  i,a,j,apla,alpha,facint);
	  #endif
	}
	*/

	Cn = 0.5*mpla*cos(0.5*xi)/(Omega*aux*(a*a)*apla*sin(xi));
	Cwe = mpla*aux/(Omega*e*(a*a)*apla);
	Cwi = 0.5*tan(0.5*xi)*mpla*cos(0.5*xi)/(Omega*aux*(a*a)*apla);
	Ce = -aux*mpla/(Omega*(a*a)*e*apla);
	Ciw = -tan(0.5*xi)*mpla/(Omega*(a*a)*aux*apla);
	Cin = -mpla/(Omega*(a*a)*aux*sin(xi)*apla);

	nodo_pla = OmM;
	w_pla = wM;
	nodo = Om;

	/*      
		c*******************************************************************      
		c    dR_D/de, dR_D/dw, dR_D/ds, dR_D/dOm:
		c*******************************************************************
		c    Terminos que NO dependen de cada resonancia:
	*/
	dc0 = 2.0*si*f3(alpha)+2.0*si*(epla*epla+e*e)*f7(alpha)+
	  4.0*(si*si*si)*f8(alpha)+2.0*si*(sipla*sipla)*f9(alpha);
	c0 = 2.0*e*f2(alpha)+2.0*e*(epla*epla)*f5(alpha)+
	  4.0*(e*e*e)*f4(alpha)+2.0*e*(sipla*sipla+si*si)*f7(alpha); 
	wp0 = Cwe*c0+Cwi*dc0;
	Op0 = Cn*dc0;

	/*
	  c*******************************************************************      
	  c    Terminos que dependen de cada resonancia:
	  c*******************************************************************      
	*/
	c1 = epla*f10(alpha) + 3.*(e*e)*epla*f11(alpha) +
	  (epla*epla*epla)*f12(alpha) +
	  epla*(si*si + sipla*sipla)*f13(alpha);
    
	dc1 = 2.0*si*e*epla*f13(alpha);

	wp[1]  = Cwe*c1 + Cwi*dc1;
	Rcos[1] = cos( w_pla - w );

	a1 = e*epla*f10(alpha) + (e*e*e)*epla*f11(alpha) +
	  (epla*epla*epla)*e*f12(alpha) +
	  (si*si + sipla*sipla)*e*epla*f13(alpha);

	xep[1] = Ce*a1;
	Rsin[1] = sin( w_pla - w );
	Op[1]  = Cn*dc1;
	xip[1] = Ciw*a1;

	/*
	  c*******************************************************************      
	  c    3. nodo_pla - nodo
	  c*******************************************************************      
	*/
	c2 = 2.0*e*si*sipla*f15(alpha);
	dc2  = sipla*f14(alpha) +
	  sipla*(epla*epla + e*e)*f15(alpha) +
	  (sipla*sipla*sipla)*f16(alpha) + 3.0*(si*si)*sipla*f16(alpha);
  
	wp[2] = Cwe*c2 + Cwi*dc2;
	Rcos[2] = cos( nodo_pla - nodo ) ;
	Op[2]  = Cn*dc2;

	a2 = si*sipla*f14(alpha) + 
	  (epla*epla + e*e)*si*sipla*f15(alpha) +
	  (sipla*sipla + si*si)*si*sipla*f16(alpha);

	xip[2] = Cin*a2;
    
	/*
	  c*******************************************************************      
	  c    4. 2 ( w_pla - w )       
	  c*******************************************************************      
	*/
	c3 = 2.0*(epla*epla)*e*f17(alpha);
	wp[3]  = Cwe*c3;
	Rcos[3] = cos(2.*( w_pla - w) );
	a3 = (epla*epla)*(e*e)*f17(alpha);
	xep[3] = 2.0*Ce*a3;
	Rsin[3] = sin( 2.0*(w_pla - w) );
	xip[3] = 2.*Ciw*a3;

	/*
	  c*******************************************************************      
	  c    5. 2 ( w - nodo )
	  c*******************************************************************      
	*/
	c4 = 2.0*e*(si*si)*f18(alpha);
	dc4 = 2.0*si*(e*e)*f18(alpha);
	wp[4]  = Cwe*c4 + Cwi*dc4;
	Rcos[4] = cos( 2.*( w - nodo ) );
	a4 = (e*e)*(si*si)*f18(alpha);
	xep[4] = -2.0*Ce*a4;
	Rsin[4] = sin( 2.*( w - nodo ) );
	Op[4]  = Cn*dc4;
	xip[4] = 2.*(Cin - Ciw)*a4;

	/*
	  c*******************************************************************      
	  c    6. w + w_pla - 2*nodo
	  c*******************************************************************      
	*/
	c5 = epla*(si*si)*f19(alpha);
	dc5 = 2.0*e*epla*si*f19(alpha);
	wp[5]  = Cwe*c5 + Cwi*dc5;
	Rcos[5] = cos(w + w_pla - 2.*nodo);
	a5 = e*epla*(si*si)*f19(alpha);
	xep[5] = -Ce*a5;
	Rsin[5] = sin(w + w_pla - 2.*nodo);
	Op[5]  = Cn*dc5;
	xip[5] = (2.*Cin - Ciw)*a5;

	/*
	  c*******************************************************************      
	  c    7. 2*(w_pla - nodo)
	  c*******************************************************************      
	*/
	dc6 = 2.0*(epla*epla)*si*f20(alpha);
	wp[6]  = Cwi*dc6;
	Rcos[6] = cos(2.*(w_pla - nodo));
	Op[6]  = Cn*dc6;
	a6 = (epla*epla)*(si*si)*f20(alpha);
	xip[6] = 2.*Cin*a6;
	Rsin[6] = sin(2.*(w_pla - nodo));

	/*
	  c*******************************************************************      
	  c    8. 2*w - nodo_pla - nodo 
	  c*******************************************************************      
	*/
	c7 = 2.0*e*si*sipla*f21(alpha);
	dc7  = (e*e)*sipla*f21(alpha);
	wp[7]  = Cwe*c7 + Cwi*dc7;
	Rcos[7] = cos(2.*w - nodo_pla - nodo);
	a7 = (e*e)*si*sipla*f21(alpha);
	xep[7] = -2.*Ce*a7;
	Rsin[7] = sin(2.*w - nodo_pla - nodo);
	Op[7]  = Cn*dc7;
	xip[7] = (Cin-2.*Ciw)*a7;
      
	/*
	  c*******************************************************************      
	  c    9. w_pla - w + nodo - nodo_pla
	  c*******************************************************************      
	*/
	c8 = epla*si*sipla*f22(alpha);
	dc8 = e*epla*sipla*f22(alpha);
	wp[8]  = Cwe*c8 + Cwi*dc8 ;
	Rcos[8] = cos(w_pla - w + nodo - nodo_pla);
	a8 = e*epla*si*sipla*f22(alpha);
	xep[8] = Ce*a8;
	Rsin[8] = sin(w_pla - w + nodo - nodo_pla);
	Op[8]  = Cn*dc8;
	xip[8] = (Ciw - Cin)*a8;

	/*
	  c*******************************************************************      
	  c    10. w_pla - w - nodo + nodo_pla
	  c*******************************************************************      
	*/
	c9 = epla*si*sipla*f23(alpha);
	dc9 = e*epla*sipla*f23(alpha);
	wp[9]  = Cwe*c9 + Cwi*dc9;
	Rcos[9] = cos(w_pla - w - nodo + nodo_pla);
	a9 = e*epla*si*sipla*f23(alpha);
	xep[9] = Ce*a9;
	Rsin[9] = sin(w_pla - w - nodo + nodo_pla);
	Op[9]  = Cn*dc9;
	xip[9] = (Cin + Ciw)*a9;

	/*
	  c*******************************************************************      
	  c    11. w + w_pla - nodo - nodo_pla        
	  c*******************************************************************      
	*/
	c10 = epla*si*sipla*f24(alpha);
	dc10 = e*epla*sipla*f24(alpha);
	wp[10]  =  Cwe*c10 + Cwi*dc10;
	Rcos[10] = cos(w + w_pla - nodo - nodo_pla);
	a10 = e*epla*si*sipla*f24(alpha);
	xep[10] = -Ce*a10;
	Rsin[10] = sin(w + w_pla - nodo - nodo_pla);
	Op[10]  = Cn*dc10;
	xip[10] = (Cin - Ciw)*a10;

	/*
	  c*******************************************************************      
	  c    12. 2 w_pla - nodo_pla - nodo       
	  c*******************************************************************      
	*/
	dc11 = (epla*epla)*sipla*f25(alpha);
	wp[11] = Cwi*dc11;
	Rcos[11] = cos(2.*w_pla - nodo_pla - nodo);
	Op[11]  = Cn*dc11;
	a11 = (epla*epla)*si*sipla*f25(alpha);
	xip[11] = Cin*a11;
	Rsin[11] = sin(2.*w_pla - nodo_pla - nodo);

	/*
	  c*******************************************************************      
	  c    13. 2 w - 2*nodo_pla
	  c*******************************************************************      
	*/
	c12 = 2.0*e*(sipla*sipla)*f18(alpha);
	wp[12]  =  Cwe*c12;
	Rcos[12] = cos(2.*w - 2.*nodo_pla);
	a12 = (e*e)*(sipla*sipla)*f18(alpha);
	xep[12] = -2.*Ce*a12;
	Rsin[12] = sin(2.*w - 2.*nodo_pla);
	xip[12] = -2.0*Ciw*a12;

	/*
	  c*******************************************************************      
	  c    14. w_pla + w - 2 nodo_pla        
	  c*******************************************************************      
	*/
	c13 = epla*(sipla*sipla)*f19(alpha);
	wp[13]  = Cwe*c13;
	Rcos[13] = cos(w_pla + w - 2.*nodo_pla);
	a13 = e*epla*(sipla*sipla)*f19(alpha);
	xep[13] = -Ce*a13;
	Rsin[13] = sin(w_pla + w - 2.*nodo_pla);
	xip[13] = -Ciw*a13;

	/*
	  c*******************************************************************      
	  c    15. 2 (w_pla - nodo_pla) 
	  c*******************************************************************      
	*/

	/*
	  c*******************************************************************      
	  c    16. 2*(nodo_pla - nodo)
	  c*******************************************************************      
	*/
	dc15 = 2.0*si*(sipla*sipla)*f26(alpha);
	wp[15] = Cwi*dc15;
	Rcos[15] = cos(2.*(nodo_pla - nodo));
	Op[15]  = Cn*dc15;
	a15 = (si*si)*(sipla*sipla)*f26(alpha);
	Rsin[15] = sin(2.*(nodo_pla - nodo));
	xip[15] = 2.0*Cin*a15;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//CALCULATE INDIVIDUAL CONTRIBUTIONS
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dxdtp1=0.0;
	dxdtp2=0.0;
	dxdtp3=Op0;
	dxdtp4=wp0;
	facint=1;
	for(ikk=1;ikk<=nkk;ikk++){
	  dxdtp1+=xep[ikk]*Rsin[ikk]*facint;
	  dxdtp2+=xip[ikk]*Rsin[ikk]*facint;
	  dxdtp3+=Op[ikk]*Rcos[ikk]*facint;
	  dxdtp4+=wp[ikk]*Rcos[ikk]*facint;
	}
	//fprintf(stdout,"\tIntegration: %e %e %e %e\n",dxdtp1,dxdtp2,dxdtp3,dxdtp4);
	dxdt[1+ki]+=dxdtp1;
	dxdt[2+ki]+=dxdtp2;
	dxdt[3+ki]+=dxdtp3;
	dxdt[4+ki]+=dxdtp4;
      }
    }
    return 0;
  }

  int secular2(double *dxdt)
  {
    return 0;
  }//END SECULAR 2 
  
};
int secularFunction(double t,const double y[],double yp[],void *param)
{
  SecularEvolution *plsys=(SecularEvolution*)param;
  plsys->update(y);
  plsys->secular(yp);
  return 0;
}

int tidalAcceleration(double t,const double y[],double dydt[],params ps)
{
  SecularEvolution *plsys=(SecularEvolution*)ps;
  int Np,Nptid;
  real wterm,w220q,X220q;
  real tauTerm;
  real tauTidal,tauTriaxial,dotEtid,dota;
  real n;
  int l=2,m=2,p=0,q;
  real chitauA;
  real R,I;
  int i;
  real M,E,f,r;

  static int ncall=1;

  //Number of planets
  Np=plsys->Np;
  Nptid=plsys->Nptid;
  
  //RESET dydt
  if(Mode>0)
    memset(dydt,0,Np*NUMVARS*sizeof(dydt[0]));

  #ifdef VERBOSE
  fprintf(stdout,"Time:%e\n",t);
  fprintf(stdout,"Number of planets:%d\n",Np);
  fprintf(stdout,"Vector:");
  fprintf_vec(stdout,"%e ",y,NUMVARS*Np);
  fprintf(stdout,"Derivatives:");
  fprintf_vec(stdout,"%e ",dydt,NUMVARS*Np);
  #endif

  int ip,k,indp;
  Body b;

  if(Mode==0 || Mode==1){
    for(ip=0;ip<Nptid;ip++){
      k=NUMVARS*ip;
      indp=plsys->Iplanetstid[ip];
      b=Bodies[indp];
      n=b.n;
    
      ////////////////////////////////
      //d theta / dt = omega
      ////////////////////////////////
      dydt[5+k]=y[6+k];
  
      ////////////////////////////////
      //Torques
      ////////////////////////////////
      //Tidal torque
      tauTidal=0;
      dotEtid=0;
      wterm=(l-2*p)*n-m*y[6+k];
      for(q=-2;q<=7;q++){
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
      M=(real)fmod((real)(n*t),(real)(2*PI));
      E=solveKepler(b.e,M);
      f=2*atan(sqrt((1+b.e)/(1-b.e))*sin(E/2)/cos(E/2));
      r=b.a*(1-b.e*cos(E));

      //Triaxial torque
      tauTriaxial=1.5*(b.B-b.A)*b.mu*sin(2*(f-y[5+k]))/(r*r*r);

      //Store values for external use
      TauTidal=tauTidal/b.C;
      TauTriax=tauTriaxial/b.C;

      ////////////////////////////////
      //d omega / dt = Torque
      ////////////////////////////////
      dydt[6+k]=(tauTidal+tauTriaxial)/b.C;

      ////////////////////////////////
      //d E / dt:
      ////////////////////////////////
      dydt[7+k]=dotEtid;

      ////////////////////////////////
      //Orbit evolution:
      ////////////////////////////////
      dydt[0+k]=0.0;
      dydt[1+k]=0.0;

      ////////////////////////////////
      //d a / dt:
      ////////////////////////////////
      dota=2*b.a*b.a/(b.mu*b.M)*dotEtid;
      dydt[0+k]=dota;

      ////////////////////////////////
      //d e / dt:
      ////////////////////////////////
      dydt[1+k]=
	1/(2*b.a)*((1-b.e*b.e)*dota/2-(1/b.M)*sqrt(b.a*(1-b.e*b.e)/b.mu)*tauTidal);

      #ifdef VERBOSE
      fprintf(stdout,"Tidal changes:\n");
      fprintf_vec(stdout,"%e ",dydt,NUMVARS*Np);
      Verbose("KEPLER: M=%e,E=%e,f=%e,r=%e\n",M,E,f,r);
      Verbose("tauTidal=%e,tauTriax=%e,dwdt=%e,dEdt=%e\n",
	      tauTidal,tauTriaxial,dydt[6+k],dydt[7+k]);
      #endif
    }
  }
  ////////////////////////////////
  //Secular Contribution
  ////////////////////////////////
  if(Mode>0){
    plsys->update(y);
    plsys->secular(dydt);
  } 

  #ifdef VERBOSE
  fprintf(stdout,"Secular changes:\n");
  fprintf_vec(stdout,"%e ",dydt,NUMVARS*Np);
  //exit(0);
  #endif

  return GSL_SUCCESS;
}
