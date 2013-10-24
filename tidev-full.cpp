#include <tidev.cpp>
#define ODEMETHOD gsl_odeiv2_step_rk4

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
  configValue(real,ftmin,"ftmin");
  configValue(real,ftmax,"ftmax");
  configValue(real,dtstep,"dtstep");
  configValue(real,dtscreen,"dtscreen");
  configValue(real,tend,"tend");
  configValue(real,thetaini,"thetaini");
  configValue(real,Pini,"Pini");

  ////////////////////////////////////////////////////////////////////////
  //ECCENTRICITY FUNCTIONS
  ////////////////////////////////////////////////////////////////////////
  readGEcc(hansen_coefs);

  ////////////////////////////////////////////////////////////////////////
  //PREPARE INTEGRATOR
  ////////////////////////////////////////////////////////////////////////
  //BODY
  IBody=1;

  //Variables (5): theta,omega,E,a,e
  gsl_odeiv2_system sys={tidalAcceleration,NULL,5,NULL};
  gsl_odeiv2_driver* driver=
    gsl_odeiv2_driver_alloc_y_new(&sys,ODEMETHOD,
				  1E-3,0,1E-6);
  
  ////////////////////////////////////////////////////////////////////////
  //RUN INTEGRATION
  ////////////////////////////////////////////////////////////////////////
  int status;
  real y[5],dydt[5];
  real accel,Etid;
  real t=0,tstep=dtstep;
  int NT=0;
  double TINI,T1,T2,TIME,TEND,TAVG;
  double Prot,dtmin,dtmax;

  //Initial conditions
  y[0]=thetaini*PI/180;	        //Angulo inicial (probabilidad de captura: P(y(1)))
  y[1]=2*PI/(Pini*HOURS/UT); 	//Initial angular velocity
  y[2]=0.0; 			//Disipated energy
  y[3]=Bodies[IBody].a;		//Semimajor axis init
  y[4]=Bodies[IBody].e;		//Eccentricity init

  //Integrate!
  file fl=fopen(outputFile,"w");
  TINI=T1=Time();
  for(tstep=dtstep;tstep<tend;tstep+=dtstep){

    //STEP SIZE
    Prot=2*PI/y[1];
    dtmin=ftmin*Prot;
    dtmax=ftmax*Prot;
    gsl_odeiv2_driver_set_hmin(driver,dtmin);
    gsl_odeiv2_driver_set_hmax(driver,dtmax);
    fprintf(stdout,"dtmin = %e, dtmax = %e\n",dtmin,dtmax);

    //STEP
    status=gsl_odeiv2_driver_apply(driver,&t,tstep,y);

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
      printf("t=%e, theta=%e, Omega/n=%e, Etid=%e, a=%e, e=%e, cpu time = %e sec, cumulative = %e sec\n",
	     t,y[0],y[1]/Bodies[IBody].n,Etid,y[3],y[4],TIME,(T2-TINI)*MICRO);
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
