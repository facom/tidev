#include <tidev.cpp>
#define ODEMETHOD gsl_odeiv2_step_rk2

int main(int argc,char *argv[])
{
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
  //INITIAL CONDITIONS PLANET AND PERTURBERS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int npert;

  //INITIAL STATE VECTOR
  double* x=(double*)calloc(5,sizeof(double));

  //GRADIENT VECTOR
  double* dxdt=(double*)calloc(5,sizeof(double));

  //PERTURBER STATE VECTOR
  double** xp=matrixAlloc(NPlanets,6);

  //LAPLACE COEFFICIENTS
  LaplaceCoefficients* laps=(LaplaceCoefficients*)calloc(NPlanets,
							 sizeof(LaplaceCoefficients));
  
  //INITIALIZATION STATE VECTOR PERTURBED BODY
  IBody=2;
  x[0]=Bodies[IBody].a;
  x[1]=Bodies[IBody].e;
  x[2]=Bodies[IBody].xi*D2R;
  x[3]=Bodies[IBody].Om*D2R;
  x[4]=fmod((Bodies[IBody].w+Bodies[IBody].Om)*D2R,2*PI);
  fprintf(stdout,"Perturbed planet %d: ",IBody);
  fprintf_vec(stdout,"%e ",x,5);

  //INITIALIZATION STATE VECTOR PERTURBER BODIES
  int j=0;
  double alpha;
  //1-b,2-c,3-d,4-e,5-f,6-g
  for(int i=1;i<=NPlanets;i++){
    if(i!=IBody && (i==4)){
      xp[j][0]=Bodies[i].a;
      xp[j][1]=Bodies[i].e;
      xp[j][2]=Bodies[i].xi*D2R;
      xp[j][3]=Bodies[i].Om*D2R;
      xp[j][4]=fmod((Bodies[i].w+Bodies[i].Om)*D2R,2*PI);
      xp[j][5]=Bodies[i].M;
      alpha=Bodies[IBody].a<xp[j][0]?Bodies[IBody].a/xp[j][0]:xp[j][0]/Bodies[IBody].a;
      laps[j].set(alpha);
      laps[j].coefcs();

      #ifdef VERBOSE
      Verbose("Perturber %d -> Planet %d:%s\n",j,i,STR(Bodies[i].name));
      Verbose("Elements:");
      for(int k=0;k<=5;k++)
	fprintf(stdout,"%e ",xp[j][k]);
      fprintf(stdout,"\n");
      #endif

      j++;
    }
  }
  npert=j;
  fprintf(stdout,"Number of perturbing planets: %d\n",npert);
  fprintf(stdout,"Mu = %e\n",GPROG*Bodies[0].M);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIALIZE SECULAR EVOLUTION SYSTEM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SecularEvolution plsys;
  plsys.set(npert,x,xp,laps);
  plsys.secular(dxdt);
  fprintf_vec(stdout,"%e ",dxdt,5);

  //exit(0);
  
  gsl_odeiv2_system sys={secularFunction,NULL,5,&plsys};

  double NPeriod=1E5;
  double P=Bodies[IBody].P;
  double t=0.0;
  double dt=P/100.0;
  double tend=NPeriod*P;
  tend=10000;
  double xout[5];

  fprintf(stdout,"Integration information:\n");
  fprintf(stdout,"\tP = %e\n",P);
  fprintf(stdout,"\tdt = %e\n",dt);
  fprintf(stdout,"\ttend = %e\n",tend);

  gsl_odeiv2_driver* d=
    gsl_odeiv2_driver_alloc_y_new(&sys,
				  ODEMETHOD,
				  dt,
				  0.0,1E-6);

  FILE *fl=fopen("results.dat","w");

  int i=0;
  do{
    gsl_odeiv2_driver_apply_fixed_step(d,&t,dt,1,x);

    if((i%((int)NPeriod))==0){
      fprintf(stdout,"t = %e\n",t);

      xout[0]=x[0];
      xout[1]=x[1];

      xout[2]=fmod(x[2]/D2R,360.0);
      if(xout[2]<0) xout[2]+=360;

      xout[3]=fmod(x[3]/D2R,360.0);
      if(xout[3]<0) xout[3]+=360;

      xout[4]=fmod(x[4]/D2R,360.0);
      xout[4]=fmod(xout[4]-xout[3],360);
      if(xout[4]<0) xout[4]+=360;

      fprintf(fl,"%e %e %e %e %e %e\n",t,xout[0],xout[1],xout[2],xout[3],xout[4]);
    }
    i++;
  }while(t<tend);

  fclose(fl);

  return 0;
}
