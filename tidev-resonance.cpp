#include <tidev.cpp>

//#define METHOD gsl_odeiv2_step_rk8pd
//#define METHOD gsl_odeiv2_step_rk2
#define METHOD gsl_odeiv2_step_rk4


#define OPTIONS "c:o:"
const char* Usage=
  MULTI(
	Usage:\n
	./tidev-resonance.cpp [-c <configuration-file>] [-o <output-file>]\n
	\n
	<configuration-file>: file with parameters.  Default: tidev.cfg\n
	<output-file>: file for results\n
	);

int main(int argc,char *argv[])
{
  ////////////////////////////////////////////////////////////////////////
  //COMMAND LINE OPTIONS
  ////////////////////////////////////////////////////////////////////////
  const char* configFile="tidev.cfg";
  const char* outputFile="output.dat";
  while((Option=getopt(argc,argv,OPTIONS))!=-1){
    switch(Option){
    case 'c':
      configFile=optarg;
      break;
    case 'o':
      outputFile=optarg;
      break;
    case '?':
      fprintf(stderr,"%s",Usage);
      exit(1);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  //INITIALIZATION
  ////////////////////////////////////////////////////////////////////////
  //TURN OFF GSL ERROR 
  errorOff();
  
  //LOAD CONFIGURATION
  configInit();
  configLoad(configFile);

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
  //BEHAVIOR
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(string,hansen_coefs,"hansen_coefs");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PHYSICAL SYSTEM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  readBodies();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //NUMERICAL
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(real,dt,"dt");
  configValue(real,dtstep,"dtstep");
  configValue(real,dtscreen,"dtscreen");
  configValue(real,tend,"tend");
  configValue(real,thetaini,"thetaini");
  configValue(real,omegaini,"omegaini");

  ////////////////////////////////////////////////////////////////////////
  //ECCENTRICITY FUNCTIONS
  ////////////////////////////////////////////////////////////////////////
  readGEcc(hansen_coefs);

  ////////////////////////////////////////////////////////////////////////
  //PREPARE INTEGRATOR
  ////////////////////////////////////////////////////////////////////////
  //Variables (5): theta,omega,E,a,e
  gsl_odeiv2_system sys={tidalAcceleration,NULL,5,NULL};
  gsl_odeiv2_driver* driver=
    gsl_odeiv2_driver_alloc_y_new(&sys,METHOD,
				  dt,0,1E-5);
  
  ////////////////////////////////////////////////////////////////////////
  //RUN INTEGRATION
  ////////////////////////////////////////////////////////////////////////
  real y[5],dydt[5];
  real accel,Etid;
  real t=0,tstep=dtstep;
  int NT=0;
  double TINI,T1,T2,TIME,TEND,TAVG;

  //Initial conditions
  y[0]=thetaini*PI/180;	        //Angulo inicial (probabilidad de captura: P(y(1)))
  y[1]=omegaini*Bodies[IBody].n; 	//Periodo inicial equivalente a 10 h 
  y[2]=0.0; 			//Disipated energy
  y[3]=Bodies[IBody].a;		//Semimajor axis init
  y[4]=Bodies[IBody].e;		//Eccentricity init

  //Integrate!
  file fl=fopen(outputFile,"w");
  TINI=T1=Time();
  for(tstep=0;tstep<tend;tstep+=dtstep){

    //STEP
    int status=gsl_odeiv2_driver_apply(driver,&t,tstep,y);

    //COMPUTE ACCELERATION AND DISSIPATED ENERGY
    tidalAcceleration(tstep,y,dydt,NULL);
    accel=dydt[1];
    Etid=y[2]*(UM*UL*UL/(UT*UT));

    //PRINT SCREEN REPORT
    fprintf(fl,"%e %e %e %e %e %e %e %e %e\n",
	    tstep,
	    y[0],y[1]/Bodies[IBody].n,Etid,y[3],y[4],
	    accel,TauTidal,TauTriax);

    //STORE IN OUTPUT FILE
    if(fmod(tstep,dtscreen)<1E-5){
      T2=Time();
      TIME=(T2-T1)*MICRO;
      printf("t=%e, y[0]=%e, y[1]/n=%e, y[2]=%e, y[3]=%e, y[4]=%e, time = %e sec\n",
	     t,y[0],y[1]/Bodies[IBody].n,Etid,y[3],y[4],TIME);
      if(t>dtstep){
	TAVG+=TIME;
	NT++;
      }
      T1=Time();
    }

  }
  TEND=Time();
  TAVG/=NT;
  printf("Average time: %e secs\n",TAVG);
  printf("Total time: %e secs\n",(TEND-TINI)*MICRO);
  fclose(fl);

  return 0;
}
