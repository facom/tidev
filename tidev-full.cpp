//NUMBER OF VARIABLES
#include <tidev.cpp>

int main(int argc,char *argv[])
{
  FILE *ff;
  int qdynamic;
  int i,j,n;
  int ip,jp;
  int k=0;
  double tini;
  char fmode[]="w";

  ////////////////////////////////////////////////////////////////////////
  //INITIALIZATION
  ////////////////////////////////////////////////////////////////////////
  //TURN OFF GSL ERROR 
  errorOff();

  //LOAD CONFIGURATION
  configLoad("tidev.cfg");
  configInit();

  //RANDOM NUMBER GENERATION
  RanGen=gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(RanGen,time(NULL));

  ////////////////////////////////////////////////////////////////////////
  //CONFIGURATION PARAMETERS
  ////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //MODE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(string,mode,"mode");
  if(!mode.compare("tidal")) Mode=0;
  else if(!mode.compare("tidal+semisecular")) Mode=1;
  else if(!mode.compare("tidal+secular")) Mode=2;
  else if(!mode.compare("secular")) Mode=3;
  else{
    fprintf(stderr,"Mode not recognized\n");
    exit(1);
  }
  char secstr[100];
  sprintf(secstr,"tidal");
  if(Mode==1)
    sprintf(secstr,"tidal+semisecular");
  if(Mode==2)
    sprintf(secstr,"tidal+secular");
  if(Mode==3)
    sprintf(secstr,"secular");
  fprintf(stdout,"Mode: %s (%d)\n",STR(mode),Mode);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //UNITS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(Real,ul,"ul");
  configValue(Real,um,"um");
  configValue(Real,ut,"ut");
  configValue(Real,uG,"uG");
  setUnits(ul,um,ut,uG);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //HANSEN COEFFICIENTS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(string,hansen_coefs,"hansen_coefs");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PHYSICAL SYSTEM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  readBodies();
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //NUMERICAL
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  configValue(Real,ftint,"ftint");
  configValue(Real,ftorb,"ftorb");
  configValue(Real,dtstep,"dtstep");
  configValue(Real,dtscreen,"dtscreen");
  configValue(Real,dtdump,"dtdump");
  configValue(Real,tint,"tint");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //LIST OF PLANETS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int iplanets[100],iplanetstid[100];
  int niplanets=0;
  int niplanetstid=0;
  fprintf(stdout,"Planets:\n");
  for(i=1;i<NBodies;i++){
    fprintf(stdout,"\tPlanet %d %s:",i,STR(Bodies[i].name));
    if(Bodies[i].active){
      fprintf(stdout,"Active ");
      iplanets[niplanets++]=i;
      if(Bodies[i].tidal && Mode<=1){
	iplanetstid[niplanetstid++]=i;
	fprintf(stdout,"Tidal");
      }
    }else{
      fprintf(stdout,"Inactive");
    }
    fprintf(stdout,"\n");
  }
  if(!niplanets){
    fprintf(stderr,"No planet active.  Please active at least one planet.\n");
    exit(1);
  }
  if(!niplanetstid && Mode<=1){
    fprintf(stderr,"No planet tidally active.  Please active at least one planet.\n");
    exit(1);
  }
  fprintf(stdout,"Active: %d, Tidal: %d\n",niplanets,niplanetstid);

  ////////////////////////////////////////////////////////////////////////
  //LOAD ECCENTRICITY DATA (IF SEMISECULAR)
  ////////////////////////////////////////////////////////////////////////
  if(Mode==1){
    for(i=0;i<niplanets;i++){
      ip=iplanets[i];
      loadEccData(&Bodies[ip]);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  //ECCENTRICITY FUNCTIONS
  ////////////////////////////////////////////////////////////////////////
  readGEcc(hansen_coefs);

  ////////////////////////////////////////////////////////////////////////
  //PREPARE INTEGRATOR
  ////////////////////////////////////////////////////////////////////////
  double* x=(double*)calloc(NUMVARS*niplanets,sizeof(double));
  double* dx=(double*)calloc(NUMVARS*niplanets,sizeof(double));
  double* xout=(double*)calloc(NUMVARS*niplanets,sizeof(double));
  double* dxdt=(double*)calloc(NUMVARS*niplanets,sizeof(double));

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIAL CONDITIONS
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  double phase;
  FILE* ft;
  if((ft=fopen("evolution.dump","r"))!=NULL){
    sprintf(fmode,"a");
    fscanf(ft,"%lf",&tini);
    fprintf(stdout,"Initial conditions from dump with tini = %e yrs.\n",tini);
    for(i=0;i<niplanets*NUMVARS;i++) fscanf(ft,"%lf",&x[i]);
    fclose(ft);
  }else{
    fprintf(stdout,"No dump.\n");
    for(i=0;i<niplanets;i++){
      k=8*i;
      ip=iplanets[i];
      //0:a
      x[0+k]=Bodies[ip].a;
      //1:e
      x[1+k]=Bodies[ip].e;
      if(Mode==1){
	//INITIAL PHASE
	/*
	phase=2*PI*Random();
	for(n=0;n<Bodies[ip].NEccFt;n++)
	  Bodies[ip].EccFt[n][5]=phase;
	//*/
	x[1+k]=fourierSeries(0.0,Bodies[ip].EccFt,Bodies[ip].NEccFt,Bodies[ip].TEccFt);
      }
      //2:Inclination
      x[2+k]=Bodies[ip].I*D2R;
      //3:Omega
      x[3+k]=Bodies[ip].Om*D2R;
      //4:w + Omega
      x[4+k]=fmod((Bodies[ip].w+Bodies[ip].Om)*D2R,2*PI);
      //5:theta
      x[5+k]=Bodies[ip].thetaini;
      //6:omega_rot
      x[6+k]=2*PI/Bodies[ip].Pini;
      //7:Etid
      x[7+k]=0.0;
    }
  }

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //CREATING OUTPUT FILES
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  char fname[100];
  FILE** fls;
  fls=(FILE**)calloc(niplanets,sizeof(FILE*));
  fprintf(stdout,"Creating output files in mode '%s':\n",fmode);
  for(i=0;i<niplanets;i++){
    //sprintf(fname,"evolution_%s-%s.dat",secstr,STR(Bodies[iplanets[i]].name));
    sprintf(fname,"evolution-%s.dat",STR(Bodies[iplanets[i]].name));
    fprintf(stdout,"\tFile: %s\n",fname);
    fls[i]=fopen(fname,fmode);
    fprintf(fls[i],"%-23s %-23s %-23s %-23s %-23s %-23s %-23s %-23s %-23s %-23s\n",
	    "#1:t","2:a","3:e","4:I","5:Om","6:w+Om","7:Theta","8:Omega/n","9:Etid","10:accel");
  }

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIALIZING LAPLACE COEFFICIENTS
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  double alpha;
  LC** Laps=laplaceAlloc(niplanets,niplanets);
  for(i=0;i<niplanets;i++){
    ip=iplanets[i];
    for(j=i+1;j<niplanets;j++){
      jp=iplanets[j];
      alpha=ALPHA(Bodies[ip].a,Bodies[jp].a);

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
  plsys.set(niplanets,iplanets,niplanetstid,iplanetstid,x,Laps);

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIALIZE SYSTEM
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  gsl_odeiv2_system sys={tidalAcceleration,NULL,NUMVARS*niplanets,&plsys};
  gsl_odeiv2_step* step=gsl_odeiv2_step_alloc(ODEMETHOD,NUMVARS*niplanets);

  ////////////////////////////////////////////////////////////////////////
  //RUN INTEGRATION
  ////////////////////////////////////////////////////////////////////////
  double t,tstep;

  int status;
  Real accel,Etid;

  double nP,Prot,Protmin,aorb,norb,Porb,Porbmin,dtint=0;
  double h;
  int nstep;

  int NT=0;
  double TINI,T1,T2,TIME,TEND,TAVG=0;
  TINI=T1=Time();

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //INITIAL REPORT
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  if(ftint==0 && ftorb==0) dtint=dtstep;
  fprintf(stdout,"Integration properties:\n");
  fprintf(stdout,"\ttini = %.2e yrs\n",tini);
  fprintf(stdout,"\t\tdtint/Pmin = %.2e\n",ftint); 
  fprintf(stdout,"\t\tdtint/Porb = %.2e\n",ftorb);
  fprintf(stdout,"\t\tdtint = %.2e yrs\n",dtint);
  fprintf(stdout,"\t\tdtstep = %.2e yrs\n",dtstep);
  fprintf(stdout,"\t\tdtscreen = %.2e yrs\n",dtscreen);
  fprintf(stdout,"\t\tdtdump = %.2e yrs\n",dtdump);
  fprintf(stdout,"\ttint = %.2e yrs\n",tint);
  fprintf(stdout,"\ttend = %.2e yrs\n",tini+tint);
  
  fprintf(stdout,"Integration starts at: tstep = %e yrs\n",tini);

  tstep=tini;

  do{
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //CHOOSE STEP SIZE
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if(ftint>0){
      Protmin=1E100;
      for(i=0;i<niplanets;i++){
	ip=iplanets[i];
	k=NUMVARS*i;
	Prot=2*PI/x[6+k];
	if((Bodies[ip].tidal && Mode<=1)||Mode==2)
	  Protmin=Prot<=Protmin?Prot:Protmin;
      }
      dtint=ftint*Protmin;
    }
    else if(ftorb>0){
      Porbmin=1E100;
      for(i=0;i<niplanets;i++){
	ip=iplanets[i];
	k=NUMVARS*i;
	aorb=x[0+k];
	Porb=Bodies[ip].P;
	Porbmin=Porb<=Porbmin?Porb:Porbmin;
      }
      dtint=ftorb*Porbmin;
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //STEP
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    nstep=(int)(dtstep/dtint);
    h=dtstep/nstep;
    t=tstep;
    do{
      status=gsl_odeiv2_step_apply(step,t,h,x,dx,NULL,dxdt,&sys);
      t+=h;
    }while((t-tstep)<dtstep);
    
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INCREASING
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    tstep+=dtstep;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //DUMP
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if(fmod(tstep,dtdump)<1E-5){
      ft=fopen("evolution.dump","w");
      fprintf(ft,"%-23.17e ",tstep);
      for(i=0;i<niplanets*NUMVARS;i++)
	fprintf(ft,"%-23.17e ",x[i]);
      for(i=0;i<niplanets*NUMVARS;i++)
	fprintf(ft,"%-23.17e ",dxdt[i]);
      fclose(ft);
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //STORE RESULT
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    for(i=0;i<niplanets;i++){
      k=NUMVARS*i;
      ip=iplanets[i];

      //Energy dissipated
      Etid=x[7+k]*(UM*UL*UL/(UT*UT));

      //Total tidal acceleration
      accel=dxdt[6+k];
      
      //Period and n
      aorb=x[0+k];
      Porb=2*PI*sqrt((aorb*aorb*aorb)/Bodies[ip].mu);
      norb=2*PI/Porb;
      Bodies[ip].P=Porb;
      Bodies[ip].n=norb;

      //Columns: 0:t,1:a(0),2:e(1),3:theta(5),4:Omega(6)/n,5:Etid(7)
      fprintf(fls[i],"%-23.17e %-23.17e %-23.17e %-23.17e %-23.17e %-23.17e %-23.17e %-23.17e %-23.17e %-23.17e\n",
	      tstep,
	      x[0+k],/*a*/
	      x[1+k],/*e*/
	      x[2+k],/*I*/
	      x[3+k],/*Om*/
	      x[4+k],/*w+Om*/
	      x[5+k],/*theta*/
	      x[6+k]/Bodies[ip].n,/*Omega/n*/
	      Etid,/*Tidal energy*/
	      accel/*Tidal total acceleration*/
	      );
      fflush(fls[i]);
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //SCREEN
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if(fmod((tstep-tini),dtscreen)<1E-5){
      T2=Time();
      TIME=(T2-T1)*MICRO;
      fprintf(stdout,
	      "\tt = %.2e, h = %.2e, nsteps = %d. Timing: Timestep = %.2e secs, Cumulative = %.2e secs\n",
	      tstep,dtint,nstep,TIME,(T2-TINI)*MICRO);
      if(tstep>(tini+dtstep)){
	TAVG+=TIME;
	NT++;
      }
      T1=Time();

      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //ACCELERATE
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      configLoad("dynamic.cfg");
      qdynamic=CFG.lookup("dynamic");
      if(qdynamic){
	ftint=CFG.lookup("ftint");
	dtstep=CFG.lookup("dtstep");
	dtscreen=CFG.lookup("dtscreen");
	dtdump=CFG.lookup("dtdump");
	tint=CFG.lookup("tint");
      }
    }
  }while((tstep-tini)<tint);
  fprintf(stdout,"Integration ends at: tstep = %e yrs\n",tstep);

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //DUMPING
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ft=fopen("evolution.dump","w");
  fprintf(ft,"%-23.17e ",tstep);
  for(i=0;i<niplanets*NUMVARS;i++) fprintf(ft,"%-23.17e ",x[i]);
  fclose(ft);

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //TIMING
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  TEND=Time();
  TAVG/=NT;
  printf("Average time: %e secs\n",TAVG);
  printf("Total time: %e secs\n",(TEND-TINI)*MICRO);
  for(i=0;i<niplanets;i++) fclose(fls[i]);

  return 0;
}
