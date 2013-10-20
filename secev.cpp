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

  ////////////////////////////////////////////////////////////////////////
  //PROGRAM
  ////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //FAKE SYSTEM
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //*
  Bodies[0].M=1.0;
  GPROG=4*PI2;

  Bodies[1].name="perturber";
  Bodies[1].M=0.33;
  Bodies[1].a=1000.0;
  Bodies[1].e=0.3;
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

  NBodies=3;
  NPlanets=2;
  //*/

  fprintf(stdout,"Number of bodies:%d\n",NBodies);
  fprintf(stdout,"Number of planets:%d\n",NPlanets);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIAL CONDITIONS PLANET AND PERTURBERS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int npert;
  double* x=(double*)calloc(5,sizeof(double));
  double* dxdt=(double*)calloc(5,sizeof(double));
  double** xp=matrixAlloc(NPlanets,6);
  LaplaceCoefficients* laps=(LaplaceCoefficients*)calloc(NPlanets,
							 sizeof(LaplaceCoefficients));
  int j=0;

  IBody=2;
  x[0]=Bodies[IBody].a;
  x[1]=Bodies[IBody].e;
  x[2]=Bodies[IBody].xi*D2R;
  x[3]=Bodies[IBody].Om*D2R;
  x[4]=fmod((Bodies[IBody].w+Bodies[IBody].Om)*D2R,2*PI);
  fprintf(stdout,"Perturbed planet %d: ",IBody);
  fprintf_vec(stdout,"%e ",x,5);

  for(int i=1;i<=NPlanets;i++){
    if(i!=IBody){
      if(Bodies[i].a<Bodies[IBody].a) continue;
      xp[j][0]=Bodies[i].a;
      xp[j][1]=Bodies[i].e;
      xp[j][2]=Bodies[i].xi*D2R;
      xp[j][3]=Bodies[i].Om*D2R;
      xp[j][4]=fmod((Bodies[i].w+Bodies[i].Om)*D2R,2*PI);
      xp[j][5]=Bodies[i].M;
      fprintf(stdout,"Perturber %d -> Planet %d:%s\n",j,i,STR(Bodies[i].name));
      fprintf(stdout,"Elements:");
      for(int k=0;k<=5;k++)
	fprintf(stdout,"%e ",xp[j][k]);
      fprintf(stdout,"\n");
      laps[j].set(Bodies[IBody].a/xp[j][0]);
      laps[j].coefcs();
      j++;
    }
  }
  npert=j;
  fprintf(stdout,"Number of perturbing planets: %d\n",npert);

  SecularEvolution plsys;
  plsys.set(npert,x,xp,laps);
  plsys.secular(dxdt);
  fprintf_vec(stdout,"%e ",dxdt,5);

  exit(0);

  double tend = 1.0E8;
  double Tp = pow(x[0],1.5);
  double dt = Tp/100.0;
  double t = 0.0;

  int jj;
  double te;
  double temp[6],t1;

  gsl_odeiv2_system sys={secularFunction,NULL,6,&plsys};
  gsl_odeiv2_driver* d=
    gsl_odeiv2_driver_alloc_y_new(&sys,
				  ODEMETHOD,
				  dt,
				  0.0,1E-6);
  FILE *fl=fopen("results.dat","w");
  int il=0;
  do{
    gsl_odeiv2_driver_apply_fixed_step(d,&t,dt,1,x);
    for(int k=3;k<=5;k++){
      jj = k-2;
      temp[jj] = x[k]/D2R;
      temp[jj]= fmod(temp[jj],360.0);
      if(temp[jj]<0.0) temp[jj] = temp[jj]+360.0;
    }
    temp[3]=temp[3]-temp[2];
    t1 = fmod(temp[3],360.0);
    if(temp[3]<0.0) temp[3]=temp[3]+360.0;

    if((il%100000)==0){
      fprintf(stdout,"t = %e\n",t);
      fprintf(fl,"%e %e %e %e %e %e\n",t,x[1],x[2],temp[1],temp[2],temp[3]);
    }
    il++;
  }while(t<tend);
  fclose(fl);
  return 0;
  //*/
}
