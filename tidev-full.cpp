//NUMBER OF VARIABLES
#define NUMVARS 8 //5 elements + theta + omega + Etid
#define ODEMETHOD gsl_odeiv2_step_rk4 //Integration method
#include <tidev.cpp>

int main(int argc,char *argv[])
{
  int i,j;

  ////////////////////////////////////////////////////////////////////////
  //INITIALIZATION
  ////////////////////////////////////////////////////////////////////////
  //TURN OFF GSL ERROR 
  errorOff();

  //LOAD CONFIGURATION
  configLoad("tidev.cfg");
  configInit();

  ////////////////////////////////////////////////////////////////////////
  //CONFIGURATION PARAMETERS
  ////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //MODE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(string,mode,"mode");
  if(!mode.compare("tidal"))
    Mode=0;
  else if(!mode.compare("tidal+secular"))
    Mode=1;
  else if(!mode.compare("secular"))
    Mode=2;
  else{
    fprintf(stderr,"Mode not recognized\n");
    exit(1);
  }
  fprintf(stdout,"Mode: %s (%d)\n",STR(mode),Mode);
  
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
  #ifdef VERBOSE
  fprintf(stdout,"Number of bodies:%d\n",NBodies);
  fprintf(stdout,"Number of planets:%d\n",NPlanets);
  #endif
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //NUMERICAL
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(real,ftmin,"ftmin");
  configValue(real,ftmax,"ftmax");
  configValue(real,dtstep,"dtstep");
  configValue(real,dtscreen,"dtscreen");
  configValue(real,tend,"tend");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //LIST OF PLANETS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int iplanets[100];
  int niplanets=0;
  for(i=0;i<NBodies;i++){
    if(Bodies[i].active)
      iplanets[niplanets++]=i;
  }
  if(!niplanets){
    fprintf(stderr,"No planet active.  Please active at least one planet.\n");
    exit(1);
  }

  ////////////////////////////////////////////////////////////////////////
  //ECCENTRICITY FUNCTIONS
  ////////////////////////////////////////////////////////////////////////
  readGEcc(hansen_coefs);

  ////////////////////////////////////////////////////////////////////////
  //PREPARE OUTPUT FILES
  ////////////////////////////////////////////////////////////////////////
  char secstr[100];
  sprintf(secstr,"tidal");
  if(Mode==1)
    sprintf(secstr,"tidal+secular");
  if(Mode==2)
    sprintf(secstr,"secular");

  FILE** fls;
  fls=(FILE**)calloc(niplanets,sizeof(FILE*));
  char fname[100];
  for(i=0;i<niplanets;i++){
    sprintf(fname,"evolution_%s-%s.dat",secstr,STR(Bodies[iplanets[i]].name));
    fls[i]=fopen(fname,"w");
  }

  ////////////////////////////////////////////////////////////////////////
  //PREPARE INTEGRATOR
  ////////////////////////////////////////////////////////////////////////
  double* x=(double*)calloc(NUMVARS*niplanets,sizeof(double));
  double* xout=(double*)calloc(NUMVARS*niplanets,sizeof(double));
  double* dxdt=(double*)calloc(NUMVARS*niplanets,sizeof(double));

  int ip,jp;
  int k=0;

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIAL CONDITIONS
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  for(i=0;i<niplanets;i++){
    k=8*i;
    ip=iplanets[i];

    //0:a
    x[0+k]=Bodies[ip].a;
    //1:e
    x[1+k]=Bodies[ip].e;
    //2:xi
    x[2+k]=Bodies[ip].xi*D2R;
    //3:Omega
    x[3+k]=Bodies[ip].Om*D2R;
    //4:w + Omega
    x[4+k]=fmod((Bodies[ip].w+Bodies[ip].Om)*D2R,2*PI);
    //5:theta
    x[5+k]=Bodies[ip].thetaini;
    //6:omega
    x[6+k]=2*PI/Bodies[ip].Pini;
    //7:Etid
    x[7+k]=0.0;
    fprintf(stdout,"Planet %d (%d,%s) state vector:",i,ip,STR(Bodies[ip].name));
    fprintf_vec(stdout,"%e ",x,(i+1)*NUMVARS,k);
  }
  fprintf(stdout,"Full state vector:");
  fprintf_vec(stdout,"%e ",x,NUMVARS*niplanets);

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIALIZE KEY VARIABLES
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIALIZE SECULAR EVOLUTION
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SecularEvolution plsys;
  plsys.set(niplanets,x,iplanets,Laps);
  /*
  plsys.secular(dxdt);
  fprintf_vec(stdout,"%e ",dxdt,NUMVARS*niplanets);
  */

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIALIZE SYSTEM
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  gsl_odeiv2_system sys={tidalAcceleration,NULL,NUMVARS*niplanets,&plsys};
  gsl_odeiv2_driver* driver=gsl_odeiv2_driver_alloc_y_new(&sys,ODEMETHOD,1E-3,0,1E-6);
  
  ////////////////////////////////////////////////////////////////////////
  //RUN INTEGRATION
  ////////////////////////////////////////////////////////////////////////
  int status;
  real accel,Etid;
  real t=0;
  int NT=0;
  double TINI,T1,T2,TIME,TEND,TAVG;
  double Prot,Protmin,dtmin,dtmax;
  double TauTotal;
    
  TINI=T1=Time();
  double tstep;
  for(tstep=dtstep;tstep<tend;tstep+=dtstep){

    //STEP SIZE
    Protmin=1E100;
    for(i=0;i<niplanets;i++){
      k=8*i;
      Prot=2*PI/x[6+k];
      Protmin=Prot<=Protmin?Prot:Protmin;
    }
    dtmin=ftmin*Protmin;
    dtmax=ftmax*Protmin;
    gsl_odeiv2_driver_set_hmin(driver,dtmin);
    gsl_odeiv2_driver_set_hmax(driver,dtmax);

    fprintf(stdout,"\nIntegrating at t = %e: dtmin = %e, dtmax = %e (Mode = %s)...\n",t,dtmin,dtmax,STR(mode));

    #ifdef VERBOSE
    fprintf(stdout,"Initial state at t = %e:\n",t);
    fprintf_vec(stdout,"%e ",x,NUMVARS*niplanets);
    #endif
  
    //STEP
    status=gsl_odeiv2_driver_apply(driver,&t,tstep,x);
    
    //COMPUTE ACCELERATION AND DISSIPATED ENERGY
    tidalAcceleration(tstep,x,dxdt,&plsys);

    #ifdef VERBOSE
    fprintf(stdout,"Last state:\n");
    fprintf_vec(stdout,"%e ",x,NUMVARS*niplanets);
    fprintf(stdout,"Last derivative:\n");
    fprintf_vec(stdout,"%e ",dxdt,NUMVARS*niplanets);
    fprintf(stdout,"Time = %e secs\n",(Time()-TINI)*MICRO);
    exit(0);
    #endif

    //PRINT SCREEN REPORT
    for(i=0;i<niplanets;i++){
      k=8*i;
      ip=iplanets[i];

      //Energy dissipated
      Etid=x[2+k]*(UM*UL*UL/(UT*UT));

      //Total tidal acceleration
      accel=dxdt[6+k];
      TauTotal=0.0;

      //Columns: 0:t,1:a(0),2:e(1),3:theta(5),4:Omega(6)/n,5:Etid(7)
      fprintf(fls[i],"%e %e %e %e %e %e\n",
	      tstep,
	      x[0+k],x[1+k],
	      x[5+k],x[6+k]/Bodies[ip].n,
	      Etid,accel);
      fflush(fls[i]);
    }

    //STORE IN OUTPUT FILE
    if(fmod(tstep,dtscreen)<1E-5){
      T2=Time();
      TIME=(T2-T1)*MICRO;

      for(i=0;i<niplanets;i++){
	k=8*i;
	ip=iplanets[i];
	
	//Energy dissipated
	Etid=x[2+k]*(UM*UL*UL/(UT*UT));
	
	//Total tidal acceleration
	accel=dxdt[6+k];
	TauTotal=0.0;
	
	//Columns: 0:t,1:a(0),2:e(1),3:theta(5),4:Omega(6)/n,5:Etid(7)
	fprintf(stdout,"Planet %d (%d,%s):\n",i,ip,STR(Bodies[ip].name));
	fprintf(stdout,"\tt = %e, a = %e, e = %e, theta = %e, Omega/n = %e, Etid = %e, accel = %e, step = %e secs, cumulative = %e secs\n",
		tstep,
		x[0+k],x[1+k],
		x[5+k],x[6+k]/Bodies[ip].n,
		Etid,accel,
		TIME,(T2-TINI)*MICRO
		);
      }

      if(t>dtstep){
	TAVG+=TIME;
	NT++;
      }
      T1=Time();
    }
    //exit(0);
  }
  TEND=Time();
  TAVG/=NT;
  printf("Average time: %e secs\n",TAVG);
  printf("Total time: %e secs\n",(TEND-TINI)*MICRO);
  for(i=0;i<niplanets;i++) fclose(fls[i]);

  return 0;
}
