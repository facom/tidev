#include <tidev.cpp>
#define ODEMETHOD gsl_odeiv2_step_rk2

int main(int argc,char *argv[])
{
  int i,j;
  ////////////////////////////////////////////////////////////////////////
  //INITIALIZATION
  ////////////////////////////////////////////////////////////////////////
  //TURN OFF GSL ERROR 
  errorOff();
  
  //LOAD CONFIGURATION
  configInit();
  configLoad("tidev.cfg");
  
  ////////////////////////////////////////////////////////////////////////
  //CONFIGURATION PARAMETERS
  ////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //UNITS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(real,ul,"ul");
  configValue(real,um,"um");
  configValue(real,ut,"ut");
  configValue(real,uG,"uG");
  setUnits(ul,um,ut,uG);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PHYSICAL SYSTEM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  readBodies();

  fprintf(stdout,"Number of bodies:%d\n",NBodies);
  fprintf(stdout,"Number of planets:%d\n",NPlanets);

  ////////////////////////////////////////////////////////////////////////
  //PROGRAM
  ////////////////////////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //FAKE SYSTEM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  /*
  Bodies[0].M=1.0;
  GPROG=4*PI2;

  Bodies[1].name="perturber";
  Bodies[1].M=0.33;
  Bodies[1].a=1000.0;
  Bodies[1].e=0.3;
  Bodies[1].xi=0.001;
  Bodies[1].xi=0.001;
  Bodies[1].Om=0.001;
  Bodies[1].w=0.001;
  
  Bodies[2].name="particle";
  Bodies[2].M=0.0;
  Bodies[2].a=50.0;
  Bodies[2].e=0.001;
  Bodies[2].xi=45.0;
  Bodies[2].Om=10.0;
  Bodies[2].w=10.0;
  Bodies[2].P=pow(Bodies[2].a,1.5);

  NBodies=3;
  NPlanets=2;
  //*/
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //SELECTING INTERACTING PLANETS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //1-b,2-c,3-d,4-e,5-f,6-g
  int iplanets[]={1,2,3,4,5,6};
  int niplanets=6;

  FILE** fl;
  fl=(FILE**)calloc(niplanets,sizeof(FILE*));
  char fname[100];
  for(i=0;i<niplanets;i++){
    sprintf(fname,"results_%s.dat",STR(Bodies[iplanets[i]].name));
    fl[i]=fopen(fname,"w");
  }
  //exit(0);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIAL CONDITIONS PLANET AND PERTURBERS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIAL STATE VECTOR
  double* x=(double*)calloc(5*niplanets,sizeof(double));
  double* xout=(double*)calloc(5*niplanets,sizeof(double));
  double* dxdt=(double*)calloc(5*niplanets,sizeof(double));

  int ip,jp;
  int k=0;
  for(i=0;i<niplanets;i++){
    ip=iplanets[i];
    x[k++]=Bodies[ip].a;
    x[k++]=Bodies[ip].e;
    x[k++]=Bodies[ip].xi*D2R;
    x[k++]=Bodies[ip].Om*D2R;
    x[k++]=fmod((Bodies[ip].w+Bodies[ip].Om)*D2R,2*PI);
  }
  fprintf(stdout,"State vector:");
  fprintf_vec(stdout,"%e ",x,5*niplanets);

  //INITIALIZING LAPLACE COEFFICIENTS
  double alpha;
  LC** Laps=laplaceAlloc(niplanets,niplanets);
  for(i=0;i<niplanets;i++){
    ip=iplanets[i];
    for(j=i+1;j<niplanets;j++){
      jp=iplanets[j];
      alpha=ALPHA(Bodies[ip].a,Bodies[jp].a);
      fprintf(stdout,"i=%d (ip=%d), j=%d (jp=%d): alpha = %e\n",i,ip,j,jp,alpha);

      Laps[i][j].set(alpha,Bodies[ip].M);
      Laps[i][j].coefcs();

      Laps[j][i].set(alpha,Bodies[jp].M);
      Laps[j][i].coefcs();
    }
  }
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIALIZE SECULAR EVOLUTION SYSTEM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SecularEvolution plsys;
  plsys.set(niplanets,x,Laps);
  plsys.secular(dxdt);
  fprintf_vec(stdout,"%e ",dxdt,5*niplanets);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INTEGRATE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  gsl_odeiv2_system sys={secularFunction,NULL,5*niplanets,&plsys};

  double NPeriod=1E5;
  double P=Bodies[iplanets[0]].P;
  double t=0.0;
  double dt=P/100.0;
  double tend=NPeriod*P;
  tend=10000;

  fprintf(stdout,"Integration information:\n");
  fprintf(stdout,"\tP = %e\n",P);
  fprintf(stdout,"\tdt = %e\n",dt);
  fprintf(stdout,"\ttend = %e\n",tend);

  gsl_odeiv2_driver* d=
    gsl_odeiv2_driver_alloc_y_new(&sys,
				  ODEMETHOD,
				  dt,
				  0.0,1E-6);

  int kip;
  //exit(0);

  do{
    gsl_odeiv2_driver_apply_fixed_step(d,&t,dt,1,x);
    
    if((i%((int)NPeriod))==0){
      fprintf(stdout,"t = %e\n",t);
      for(ip=0;ip<niplanets;ip++){
	kip=5*ip;

	xout[0+kip]=x[0+kip];
	xout[1+kip]=x[1+kip];
	
	xout[2+kip]=fmod(x[2+kip]/D2R,360.0);
	if(xout[2+kip]<0) xout[2+kip]+=360;
	
	xout[3+kip]=fmod(x[3+kip]/D2R,360.0);
	if(xout[3+kip]<0) xout[3+kip]+=360;
	
	xout[4+kip]=fmod(x[4+kip]/D2R,360.0);
	xout[4+kip]=fmod(xout[4+kip]-xout[3+kip],360);
	if(xout[4+kip]<0) xout[4+kip]+=360;
	fprintf(fl[ip],"%e %e %e %e %e %e\n",t,
		xout[0+kip],xout[1+kip],xout[2+kip],xout[3+kip],xout[4+kip]);
      }
    }
    i++;
  }while(t<tend);

  for(ip=0;ip<niplanets;ip++)
    fclose(fl[ip]);

  return 0;
}
