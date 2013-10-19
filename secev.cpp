#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
using namespace std;

double df ;
double df0;
double df1, df2, df3, df4;
double df01, df02;
double df20, df22;
double df10, df11, df12;
double df13, df14;
double df31, df32, df31b, df32b;
double df100, df101, df102, df103, df104;
double aM, eM, xiM, OmM, wM, mp;

void fprintf_vec(FILE* stream,const char* fmt,const double x[],int end,int ini=0)
{
  for(int i=ini;i<end;i++){
    fprintf(stream,fmt,x[i]);
  }
  fprintf(stream,"\n");
}

double func(double th,void *param)
{
  double* ps=(double*) param;
  double j=ps[0];
  double s=ps[1];
  double alpha=ps[2];

  double f=(1/pow((1.0 - 2.0*alpha*cos(th)+alpha*alpha),s))*cos(j*th);
  
  return f;
}

gsl_integration_workspace* INTWORK;
gsl_function FUNC;
double LAPLPARAM[3];

double bs(int j,double s,double alpha)        
{
  double f,error;
  double Pi=M_PI;
  double j2=j;
  double s2=s;
  double alpha2=alpha;
  double a=0.0;
  double b=Pi;

  LAPLPARAM[0]=j;
  LAPLPARAM[1]=s;
  LAPLPARAM[2]=alpha;
  gsl_integration_qags(&FUNC,a,b,0.0,1E-6,1000,INTWORK,&s,&error);
  f=(2/M_PI)*s;
  
  return f;
}

double dbs(int j,double s,double alpha)
{
  double f,f1,f2,f3;
  
  f1=bs(j,s+1.0,alpha);
  f2=bs(j-1.0,s+1.0,alpha);
  f3=bs(j+1.0,s+1.0,alpha);
  f=s*(f2-2.0*alpha*f1+f3);

  return f;
}

double d2bs(int j,double s,double alpha)
{
  double f;
  double f0=bs(j,s+1.0,alpha);
  double f1=dbs(j,s+1.0,alpha);
  double f2=dbs(j-1.0,s+1.0,alpha);
  double f3=dbs(j+1.0,s+1.0,alpha);

  f=s*(f2-2.0*alpha*f1-2.0*f0+f3);

  return f;
}

double d3bs(int j,double s,double alpha)
{
  double f;

  double f0=dbs(j,s+1.0,alpha);
  double f1=d2bs(j,s+1.0,alpha);
  double f2=d2bs(j-1.0,s+1.0,alpha);
  double f3=d2bs(j+1.0,s+1.0,alpha);
  f=s*(f2-4.0*f0-2.0*alpha*f1+f3);

  return f;
}

double d4bs(int j,double s,double alpha)
{
  double f;

  double f0=d2bs(j,s+1.0,alpha);
  double f1=d3bs(j,s+1.0,alpha);
  double f2=d3bs(j-1.0,s+1.0,alpha);
  double f3=d3bs(j+1.0,s+1.0,alpha);
  f=s*(f2-6.0*f0-2.0*alpha*f1+f3);

  return f;
}

void coefcs(double apla,double a)
{
  double alpha=a<=apla?a/apla:apla/a;

  df=bs(0.0,0.5,alpha);

  df0=bs(1.0,1.5,alpha);
  df01=dbs(1.0,1.5,alpha);
  df02=d2bs(1.0,1.5,alpha);
  printf("df,df0,df01,df02 = %e %e %e %e\n",df,df0,df01,df02);

  df1=dbs(0.0,0.5,alpha);
  df2=d2bs(0.0,0.5,alpha);
  df3=d3bs(0.0,0.5,alpha);
  df4=d4bs(0.0,0.5,alpha);      
  printf("df1,df2,df3,df4 = %e %e %e %e\n",df1,df2,df3,df4);

  df20=bs(0.0,2.5,alpha);
  df22=bs(2.0,2.5,alpha);
  printf("df20,df22 = %e %e\n",df20,df22);

  df10=bs(1.0,0.5,alpha);
  df11=dbs(1.0,0.5,alpha); 
  df12=d2bs(1.0,0.5,alpha);
  df13=d3bs(1.0,0.5,alpha);
  df14=d4bs(1.0,0.5,alpha);
  printf("df10,df11,df12,df13,df14 = %e %e %e %e %e\n",df10,df11,df12,df13,df14);
    
  df31=dbs(0.0,1.5,alpha);
  df32=d2bs(0.0,1.5,alpha);
  df31b=dbs(2.0,1.5,alpha);
  df32b=d2bs(2.0,1.5,alpha);
  printf("df31,df32,df31b,df32b = %e %e %e %e\n",df31,df32,df31b,df32b);
    
  df100=bs(2.0,0.5,alpha);
  df101=dbs(2.0,0.5,alpha);
  df102=d2bs(2.0,0.5,alpha);
  df103=d3bs(2.0,0.5,alpha);
  df104=d4bs(2.0,0.5,alpha);
  printf("df100,df101,df102,df103,df104 = %e %e %e %e %e\n",df100,df101,df102,df103,df104);
}

double f2(double alpha){
  return alpha*df0/8.0;
}

double f3(double alpha){
  return -0.5*alpha*df0;
}

double f4(double alpha){
  return (4.0*(alpha*alpha*alpha)*df3 + (alpha*alpha*alpha*alpha)*df4)/128.0;
}

double f5(double alpha){
  return (4.0*alpha*df1 + 14.0*(alpha*alpha)*df2 +
	  8.0*(alpha*alpha*alpha)*df3 + (alpha*alpha*alpha*alpha)*df4)/32.0;
}

double f6(double alpha){
  return (24.0*alpha*df1 + 36.0*(alpha*alpha)*df2 +
	  12.0*(alpha*alpha*alpha)*df3 + (alpha*alpha*alpha*alpha)*df4)/128.0;
}

double f7(double alpha){
  return (-2.0*alpha*df0 - 4.0*(alpha*alpha)*df01 -
	  (alpha*alpha*alpha)*df02)/8.0;
}

double f8(double alpha){
  return (3.0/8.0)*(alpha*alpha)*df22 +
    (3.0/4.0)*(alpha*alpha)*df20;
}

double f9(double alpha){
  return 0.5*alpha*df0 + (3.0/4.0)*(alpha*alpha)*df22 +
    (15.0/4.0)*(alpha*alpha)*df20;
}

double f10(double alpha){
  return 0.25*(2.0*df10 - 2.0*alpha*df11 -
	       (alpha*alpha)*df12);
}

double f11(double alpha){
  return (-4.0*(alpha*alpha)*df12 - 6.0*(alpha*alpha*alpha)*df13 -
	  (alpha*alpha*alpha*alpha)*df14)/32.0;
}

double f12(double alpha){
  return (4.0*df10 - 4.0*alpha*df11 - 22.0*(alpha*alpha)*df12 -
	  10.0*(alpha*alpha*alpha)*df13 - (alpha*alpha*alpha*alpha)*df14)/32.0;
}

double f13(double alpha){
  return 0.5*(alpha*alpha)*(df31 + df31b) +
    (alpha*alpha*alpha)*(df32 + df32b)/8.0;
}

double f14(double alpha){
  return alpha*df0;
}

double f15(double alpha){
  return 0.5*alpha*df0 + (alpha*alpha)*df01 +
    0.25*(alpha*alpha*alpha)*df02;
}

double f16(double alpha){
  return -0.5*alpha*df0 - 3.0*(alpha*alpha)*df20 -
    1.5*(alpha*alpha)*df22;
}

double f17(double alpha){
  return (12.0*df100 - 12.0*alpha*df101 + 
	  6.0*(alpha*alpha)*df102 +
	  8.0*(alpha*alpha*alpha)*df103 + (alpha*alpha*alpha*alpha)*df104)/64.0;
}

double f18(double alpha){
  return 0.75*alpha*df0 + 0.5*(alpha*alpha)*df01 +
    (alpha*alpha*alpha)*df02/16.0;
}

double f19(double alpha){
  return -0.5*(alpha*alpha)*df31 - (alpha*alpha*alpha)*df32/8.0;
}

double f20(double alpha){
  return (alpha*alpha*alpha)*df02/16.0;
}

double f21(double alpha){
  return -1.5*alpha*df0 - (alpha*alpha)*df01 -
    (alpha*alpha*alpha)*df02/8.0;
}

double f22(double alpha){
  return -(alpha*alpha)*df31 - (alpha*alpha*alpha)*df32/4.0;
}

double f23(double alpha){
  return -(alpha*alpha )*df31b - (alpha*alpha*alpha)*df32b/4.0;
}

double f24(double alpha){
  return (alpha*alpha)*df31 + (alpha*alpha*alpha)*df32/4.0;
}

double f25(double alpha){
  return - (alpha*alpha*alpha)*df02/8.0;
}

double f26(double alpha){
  return 0.5*alpha*df0 + 0.75*(alpha*alpha)*df20 +
    1.5*(alpha*alpha)*df22;
}

double secular(double a,double e,double xi,double w,double Om,
	       double *wpres,double *Ompres,double *ep,double *ip)
{
  /*
c     Pert. INTERNA O EXTERNA 
c     d/dt de e, i, \varpi, \Omega 
c     4o orden en e / si = sin((1/2)*xi)           
  */
  
  double pi = M_PI;
  double deg = pi/180.0;
  double G = 4.0*(pi*pi);
  double xMc = 1.0;

  /*
      aux = dsqrt(1.0-e**2)
      si = dsin(xi/2.0)
      sini = dsin(xi)
      Omega = dsqrt(G*xMc/a**3)  
  */

  double aux = sqrt(1.0-e*e);
  double si = sin(xi/2.0);
  double sini = sin(xi);
  double Omega = sqrt(G*xMc/(a*a*a));

  /*
      wp0 = 0.0
      Op0 = 0.0
      nkk = 15
      do kk=1,nkk
        wp(kk) = 0.0
        Op(kk)= 0.0
        xep(kk)= 0.0
        xip(kk) = 0.0
        Rcos(kk) = 0.0
        Rsin(kk) = 0.0
      enddo
  */
  double wp0 = 0.0;
  double Op0 = 0.0;
  int nkk = 15;
  double wp[16],Op[16],xep[16],xip[16],Rcos[16],Rsin[16];
  for(int kk=1;kk<=nkk;kk++){
    wp[kk] = 0.0;
    Op[kk]= 0.0;
    xep[kk]= 0.0;
    xip[kk] = 0.0;
    Rcos[kk] = 0.0;
    Rsin[kk] = 0.0;
  }

  /*
      mpla = G*mp 
      apla = aM 
      epla = eM      
      sipla = dsin(xiM/2.0)
  */
  double mpla = G*mp;
  double apla = aM;
  double epla = eM;
  double sipla = sin(xiM/2.0);

  /*
      if (a.le.apla) then
        alpha = a/apla
      elseif (a.ge.apla) then
        alpha = apla/a
      endif
  */
  double alpha=a<=apla?a/apla:apla/a;

  /*
      Cn = 0.5*mpla*dcos(0.5*xi)/(Omega*aux*(a**2)*apla*dsin(xi))
      Cwe = mpla*aux/(Omega*e*(a**2)*apla)
      Cwi = 0.5*dtan(0.5*xi)*mpla*dcos(0.5*xi)/
     &      (Omega*aux*(a**2)*apla)
      Ce = -aux*mpla/(Omega*(a**2)*e*apla)
      Ciw = -dtan(0.5*xi)*mpla/(Omega*(a**2)*aux*apla)
      Cin = -mpla/(Omega*(a**2)*aux*dsin(xi)*apla)
  */
  double Cn = 0.5*mpla*cos(0.5*xi)/(Omega*aux*(a*a)*apla*sin(xi));
  double Cwe = mpla*aux/(Omega*e*(a*a)*apla);
  double Cwi = 0.5*tan(0.5*xi)*mpla*cos(0.5*xi)/(Omega*aux*(a*a)*apla);
  double Ce = -aux*mpla/(Omega*(a*a)*e*apla);
  double Ciw = -tan(0.5*xi)*mpla/(Omega*(a*a)*aux*apla);
  double Cin = -mpla/(Omega*(a*a)*aux*sin(xi)*apla);

  /*
      nodo_pla = OmM
      w_pla = wM
      nodo = Om
  */
  double nodo_pla = OmM;
  double w_pla = wM;
  double nodo = Om;

  /*      
c*******************************************************************      
c    dR_D/de, dR_D/dw, dR_D/ds, dR_D/dOm:
c*******************************************************************
c    Terminos que NO dependen de cada resonancia:
  */
  
  /*
c    1. 0       
      dc0 = 2.0*si*f3(alpha) 
     &  + 2.0*si*(epla**2 + e**2)*f7(alpha) +
     &  4.0*(si**3)*f8(alpha) + 2.0*si*(sipla**2)*f9(alpha)

      c0 = 2.0*e*f2(alpha) + 2.0*e*(epla**2)*f5(alpha) +
     &  4.0*(e**3)*f4(alpha) + 2.0*e*(sipla**2 + si**2)*f7(alpha) 

      wp0 = Cwe*c0 + Cwi*dc0 

      Op0 = Cn*dc0
  */
  double dc0 = 2.0*si*f3(alpha)+2.0*si*(epla*epla+e*e)*f7(alpha)+
    4.0*(si*si*si)*f8(alpha)+2.0*si*(sipla*sipla)*f9(alpha);
  double c0 = 2.0*e*f2(alpha)+2.0*e*(epla*epla)*f5(alpha)+
    4.0*(e*e*e)*f4(alpha)+2.0*e*(sipla*sipla+si*si)*f7(alpha); 
  wp0 = Cwe*c0+Cwi*dc0;
  Op0 = Cn*dc0;

  /*
c*******************************************************************      
c    Terminos que dependen de cada resonancia:
c*******************************************************************      
  */

  /*
c    2. w_pla - w    
      c1 = epla*f10(alpha) + 3.*(e**2)*epla*f11(alpha) +
     &  (epla**3)*f12(alpha) +
     &  epla*(si**2 + sipla**2)*f13(alpha)

      dc1 = 2.0*si*e*epla*f13(alpha)

      wp(1)  = Cwe*c1 + Cwi*dc1
      Rcos(1) = dcos( w_pla - w )

      a1 = e*epla*f10(alpha) + (e**3)*epla*f11(alpha) +
     &  (epla**3)*e*f12(alpha) +
     &  (si**2 + sipla**2)*e*epla*f13(alpha)
  */
  double c1 = epla*f10(alpha) + 3.*(e*e)*epla*f11(alpha) +
    (epla*epla*epla)*f12(alpha) +
    epla*(si*si + sipla*sipla)*f13(alpha);
    
  double dc1 = 2.0*si*e*epla*f13(alpha);

  wp[1]  = Cwe*c1 + Cwi*dc1;
  Rcos[1] = cos( w_pla - w );

  double a1 = e*epla*f10(alpha) + (e*e*e)*epla*f11(alpha) +
    (epla*epla*epla)*e*f12(alpha) +
    (si*si + sipla*sipla)*e*epla*f13(alpha);

  /*
      xep(1) = Ce*a1
      Rsin(1) = dsin( w_pla - w )      
      Op(1)  = Cn*dc1
      xip(1) = Ciw*a1
  */
  xep[1] = Ce*a1;
  Rsin[1] = sin( w_pla - w );
  Op[1]  = Cn*dc1;
  xip[1] = Ciw*a1;

  /*
c*******************************************************************      
c    3. nodo_pla - nodo
c*******************************************************************      
  */

  /*
      c2 = 2.0*e*si*sipla*f15(alpha)
      dc2  = sipla*f14(alpha) +
     &  sipla*(epla**2 + e**2)*f15(alpha) +
     &  (sipla**3)*f16(alpha) + 3.0*(si**2)*sipla*f16(alpha)

      wp(2) = Cwe*c2 + Cwi*dc2
      Rcos(2) = dcos( nodo_pla - nodo ) 

      Op(2)  = Cn*dc2

      a2 = si*sipla*f14(alpha) + 
     $ (epla**2 + e**2)*si*sipla*f15(alpha) +
     $ (sipla**2 + si**2)*si*sipla*f16(alpha) 
      xip(2) = Cin*a2
  */
  double c2 = 2.0*e*si*sipla*f15(alpha);
  double dc2  = sipla*f14(alpha) +
    sipla*(epla*epla + e*e)*f15(alpha) +
    (sipla*sipla*sipla)*f16(alpha) + 3.0*(si*si)*sipla*f16(alpha);
  
  wp[2] = Cwe*c2 + Cwi*dc2;
  Rcos[2] = cos( nodo_pla - nodo ) ;
  Op[2]  = Cn*dc2;

  double a2 = si*sipla*f14(alpha) + 
    (epla*epla + e*e)*si*sipla*f15(alpha) +
    (sipla*sipla + si*si)*si*sipla*f16(alpha);

  xip[2] = Cin*a2;
    
  /*
c*******************************************************************      
c    4. 2 ( w_pla - w )       
c*******************************************************************      
  */
  
  /*
      c3 = 2.0*(epla**2)*e*f17(alpha)
      wp(3)  = Cwe*c3
      Rcos(3) = dcos(2.*( w_pla - w) )
      a3 = (epla**2)*(e**2)*f17(alpha)
      xep(3) = 2.0*Ce*a3
      Rsin(3) = dsin( 2.0*(w_pla - w) )
      xip(3) = 2.*Ciw*a3      
  */
  double c3 = 2.0*(epla*epla)*e*f17(alpha);
  wp[3]  = Cwe*c3;
  Rcos[3] = cos(2.*( w_pla - w) );
  double a3 = (epla*epla)*(e*e)*f17(alpha);
  xep[3] = 2.0*Ce*a3;
  Rsin[3] = sin( 2.0*(w_pla - w) );
  xip[3] = 2.*Ciw*a3;

  /*
c*******************************************************************      
c    5. 2 ( w - nodo )
c*******************************************************************      
  */

  /*
      c4 = 2.0*e*(si**2)*f18(alpha)
      dc4 = 2.0*si*(e**2)*f18(alpha)
      wp(4)  = Cwe*c4 + Cwi*dc4
      Rcos(4) = cos( 2.*( w - nodo ) )
      a4 = (e**2)*(si**2)*f18(alpha)
      xep(4) = -2.0*Ce*a4
      Rsin(4) = sin( 2.*( w - nodo ) )
      Op(4)  = Cn*dc4
      xip(4) = 2.*(Cin - Ciw)*a4
  */
  double c4 = 2.0*e*(si*si)*f18(alpha);
  double dc4 = 2.0*si*(e*e)*f18(alpha);
  wp[4]  = Cwe*c4 + Cwi*dc4;
  Rcos[4] = cos( 2.*( w - nodo ) );
  double a4 = (e*e)*(si*si)*f18(alpha);
  xep[4] = -2.0*Ce*a4;
  Rsin[4] = sin( 2.*( w - nodo ) );
  Op[4]  = Cn*dc4;
  xip[4] = 2.*(Cin - Ciw)*a4;

  /*
c*******************************************************************      
c    6. w + w_pla - 2*nodo
c*******************************************************************      
  */

  /*
      c5 = epla*(si*si)*f19(alpha)
      dc5 = 2.0*e*epla*si*f19(alpha)
      wp(5)  = Cwe*c5 + Cwi*dc5
      Rcos(5) = cos(w + w_pla - 2.*nodo)
      a5 = e*epla*(si*si)*f19(alpha)
      xep(5) = -Ce*a5
      Rsin(5) = sin(w + w_pla - 2.*nodo)
      Op(5)  = Cn*dc5
      xip(5) = (2.*Cin - Ciw)*a5
  */
  double c5 = epla*(si*si)*f19(alpha);
  double dc5 = 2.0*e*epla*si*f19(alpha);
  wp[5]  = Cwe*c5 + Cwi*dc5;
  Rcos[5] = cos(w + w_pla - 2.*nodo);
  double a5 = e*epla*(si*si)*f19(alpha);
  xep[5] = -Ce*a5;
  Rsin[5] = sin(w + w_pla - 2.*nodo);
  Op[5]  = Cn*dc5;
  xip[5] = (2.*Cin - Ciw)*a5;

  /*
c*******************************************************************      
c    7. 2*(w_pla - nodo)
c*******************************************************************      
  */

  /*
      dc6 = 2.0*(epla**2)*si*f20(alpha)
      wp(6)  = Cwi*dc6
      Rcos(6) = cos(2.*(w_pla - nodo))
      Op(6)  = Cn*dc6
      a6 = (epla**2)*(si*si)*f20(alpha) 
      xip(6) = 2.*Cin*a6
	Rsin(6) = sin(2.*(w_pla - nodo))
  */
  double dc6 = 2.0*(epla*epla)*si*f20(alpha);
  wp[6]  = Cwi*dc6;
  Rcos[6] = cos(2.*(w_pla - nodo));
  Op[6]  = Cn*dc6;
  double a6 = (epla*epla)*(si*si)*f20(alpha);
  xip[6] = 2.*Cin*a6;
  Rsin[6] = sin(2.*(w_pla - nodo));

  /*
c*******************************************************************      
c    8. 2*w - nodo_pla - nodo 
c*******************************************************************      
  */
  
  /*
      c7 = 2.0*e*si*sipla*f21(alpha)
      dc7  = (e*e)*sipla*f21(alpha)
      wp(7)  = Cwe*c7 + Cwi*dc7
      Rcos(7) = cos(2.*w - nodo_pla - nodo) 
      a7 = (e*e)*si*sipla*f21(alpha)
      xep(7) = -2.*Ce*a7
      Rsin(7) = sin(2.*w - nodo_pla - nodo)
      Op(7)  = Cn*dc7
      xip(7) = (Cin-2.*Ciw)*a7
  */
  double c7 = 2.0*e*si*sipla*f21(alpha);
  double dc7  = (e*e)*sipla*f21(alpha);
  wp[7]  = Cwe*c7 + Cwi*dc7;
  Rcos[7] = cos(2.*w - nodo_pla - nodo);
  double a7 = (e*e)*si*sipla*f21(alpha);
  xep[7] = -2.*Ce*a7;
  Rsin[7] = sin(2.*w - nodo_pla - nodo);
  Op[7]  = Cn*dc7;
  xip[7] = (Cin-2.*Ciw)*a7;
      
  /*
c*******************************************************************      
c    9. w_pla - w + nodo - nodo_pla
c*******************************************************************      
  */

  /*
      c8 = epla*si*sipla*f22(alpha)
      dc8 = e*epla*sipla*f22(alpha)
      wp(8)  = Cwe*c8 + Cwi*dc8 
      Rcos(8) = cos(w_pla - w + nodo - nodo_pla)
      a8 = e*epla*si*sipla*f22(alpha)
      xep(8) = Ce*a8
      Rsin(8) = sin(w_pla - w + nodo - nodo_pla)
      Op(8)  = Cn*dc8
      xip(8) = (Ciw - Cin)*a8
  */
  double c8 = epla*si*sipla*f22(alpha);
  double dc8 = e*epla*sipla*f22(alpha);
  wp[8]  = Cwe*c8 + Cwi*dc8 ;
  Rcos[8] = cos(w_pla - w + nodo - nodo_pla);
  double a8 = e*epla*si*sipla*f22(alpha);
  xep[8] = Ce*a8;
  Rsin[8] = sin(w_pla - w + nodo - nodo_pla);
  Op[8]  = Cn*dc8;
  xip[8] = (Ciw - Cin)*a8;

  /*
c*******************************************************************      
c    10. w_pla - w - nodo + nodo_pla
c*******************************************************************      
  */

  /*
      c9 = epla*si*sipla*f23(alpha)
      dc9 = e*epla*sipla*f23(alpha)
      wp(9)  = Cwe*c9 + Cwi*dc9
      Rcos(9) = cos(w_pla - w - nodo + nodo_pla)
      a9 = e*epla*si*sipla*f23(alpha)
      xep(9) = Ce*a9
      Rsin(9) = sin(w_pla - w - nodo + nodo_pla)      
      Op(9)  = Cn*dc9
      xip(9) = (Cin + Ciw)*a9
  */
  double c9 = epla*si*sipla*f23(alpha);
  double dc9 = e*epla*sipla*f23(alpha);
  wp[9]  = Cwe*c9 + Cwi*dc9;
  Rcos[9] = cos(w_pla - w - nodo + nodo_pla);
  double a9 = e*epla*si*sipla*f23(alpha);
  xep[9] = Ce*a9;
  Rsin[9] = sin(w_pla - w - nodo + nodo_pla);
  Op[9]  = Cn*dc9;
  xip[9] = (Cin + Ciw)*a9;

  /*
c*******************************************************************      
c    11. w + w_pla - nodo - nodo_pla        
c*******************************************************************      
  */

  /*
      c10 = epla*si*sipla*f24(alpha)
      dc10 = e*epla*sipla*f24(alpha)
      wp(10)  =  Cwe*c10 + Cwi*dc10
      Rcos(10) = cos(w + w_pla - nodo - nodo_pla) 
      a10 = e*epla*si*sipla*f24(alpha)
      xep(10) = -Ce*a10
      Rsin(10) = sin(w + w_pla - nodo - nodo_pla)
      Op(10)  = Cn*dc10
      xip(10) = (Cin - Ciw)*a10
  */
  double c10 = epla*si*sipla*f24(alpha);
  double dc10 = e*epla*sipla*f24(alpha);
  wp[10]  =  Cwe*c10 + Cwi*dc10;
  Rcos[10] = cos(w + w_pla - nodo - nodo_pla);
  double a10 = e*epla*si*sipla*f24(alpha);
  xep[10] = -Ce*a10;
  Rsin[10] = sin(w + w_pla - nodo - nodo_pla);
  Op[10]  = Cn*dc10;
  xip[10] = (Cin - Ciw)*a10;

  /*
c*******************************************************************      
c    12. 2 w_pla - nodo_pla - nodo       
c*******************************************************************      
  */

  /*
      dc11 = (epla*epla)*sipla*f25(alpha)
      wp(11) = Cwi*dc11
      Rcos(11) = cos(2.*w_pla - nodo_pla - nodo)
      Op(11)  = Cn*dc11
      a11 = (epla*epla)*si*sipla*f25(alpha) 
      xip(11) = Cin*a11
	Rsin(11) = sin(2.*w_pla - nodo_pla - nodo)
  */
  double dc11 = (epla*epla)*sipla*f25(alpha);
  wp[11] = Cwi*dc11;
  Rcos[11] = cos(2.*w_pla - nodo_pla - nodo);
  Op[11]  = Cn*dc11;
  double a11 = (epla*epla)*si*sipla*f25(alpha);
  xip[11] = Cin*a11;
  Rsin[11] = sin(2.*w_pla - nodo_pla - nodo);

  /*
c*******************************************************************      
c    13. 2 w - 2*nodo_pla
c*******************************************************************      
  */

  /*
      c12 = 2.0*e*(sipla**2)*f18(alpha)
      wp(12)  =  Cwe*c12 
      Rcos(12) = cos(2.*w - 2.*nodo_pla)
      a12 = (e*e)*(sipla*sipla)*f18(alpha)
      xep(12) = -2.*Ce*a12
      Rsin(12) = sin(2.*w - 2.*nodo_pla)
      xip(12) = -2.0*Ciw*a12
  */
  double c12 = 2.0*e*(sipla*sipla)*f18(alpha);
  wp[12]  =  Cwe*c12;
  Rcos[12] = cos(2.*w - 2.*nodo_pla);
  double a12 = (e*e)*(sipla*sipla)*f18(alpha);
  xep[12] = -2.*Ce*a12;
  Rsin[12] = sin(2.*w - 2.*nodo_pla);
  xip[12] = -2.0*Ciw*a12;

  /*
c*******************************************************************      
c    14. w_pla + w - 2 nodo_pla        
c*******************************************************************      
  */

  /*
      c13 = epla*(sipla*sipla)*f19(alpha) 
      wp(13)  = Cwe*c13
      Rcos(13) = cos(w_pla + w - 2.*nodo_pla) 
      a13 = e*epla*(sipla*sipla)*f19(alpha)
      xep(13) = -Ce*a13
      Rsin(13) = sin(w_pla + w - 2.*nodo_pla)
      xip(13) = -Ciw*a13
  */
  double c13 = epla*(sipla*sipla)*f19(alpha);
  wp[13]  = Cwe*c13;
  Rcos[13] = cos(w_pla + w - 2.*nodo_pla);
  double a13 = e*epla*(sipla*sipla)*f19(alpha);
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

  /*
      dc15 = 2.0*si*(sipla*sipla)*f26(alpha)
      wp(15) = Cwi*dc15
      Rcos(15) = cos(2.*(nodo_pla - nodo))     
      Op(15)  = Cn*dc15
      a15 = (si*si)*(sipla*sipla)*f26(alpha)
      Rsin(15) = sin(2.*(nodo_pla - nodo))     
      xip(15) = 2.0*Cin*a15
  */
  double dc15 = 2.0*si*(sipla*sipla)*f26(alpha);
  wp[15] = Cwi*dc15;
  Rcos[15] = cos(2.*(nodo_pla - nodo));
  Op[15]  = Cn*dc15;
  double a15 = (si*si)*(sipla*sipla)*f26(alpha);
  Rsin[15] = sin(2.*(nodo_pla - nodo));
  xip[15] = 2.0*Cin*a15;

  /*
c*************************************************************         
      wpres = wp0
      Ompres = Op0
      ep = 0.0
      ip = 0.0
      do ikk= 1,nkk
        wpres = wpres + wp(ikk)*Rcos(ikk)
        Ompres = Ompres + Op(ikk)*Rcos(ikk)
        ep = ep + xep(ikk)*Rsin(ikk)
        ip = ip + xip(ikk)*Rsin(ikk)
      enddo
  */
  *wpres = wp0;
  *Ompres = Op0;
  *ep = 0.0;
  *ip = 0.0;
  for(int ikk=1;ikk<=nkk;ikk++){
    *wpres = *wpres + wp[ikk]*Rcos[ikk];
    *Ompres = *Ompres + Op[ikk]*Rcos[ikk];
    *ep = *ep + xep[ikk]*Rsin[ikk];
    *ip = *ip + xip[ikk]*Rsin[ikk];
  }
}

int fcn(double t,const double y[],double yp[],void *param)
{
  /*
      common /planet/ aM, eM, xiM, OmM, wM, mp  
c 
      nc = 5
      do i=1,nc
        yp(i) = 0.0d0
      enddo  
c
  */
  double a=y[1];
  double e=y[2];
  double xi=y[3];
  double Om=y[4];
  double w=y[5];
  double ep,ip,Ompres,wpres;
  
  //printf("y:");
  //fprintf_vec(stdout,"%e ",y,6);
  secular(a,e,xi,w,Om,&wpres,&Ompres,&ep,&ip);
  yp[1]=0.0;
  yp[2]=ep;
  yp[3]=ip;
  yp[4]=Ompres;
  yp[5]=wpres;
  //printf("yp:");
  //fprintf_vec(stdout,"%e ",yp,6);

  return 0;
}

int main(int argc,char *argv[])
{
  int nc = 5;
  double pi=M_PI;
  double deg = pi/180.0;
  double sec = deg/3600.0;
    
  double tfin = 1.0E8;
  
  FILE *fin;
  char line[100],filename[100];
  int n;
  double ip,a0,e0,xi0,Om0,w0;

  INTWORK=gsl_integration_workspace_alloc(1000);
  FUNC.function=&func;
  FUNC.params=&LAPLPARAM;

  //READ INPUT FILE
  fin=fopen("ersec.in","r");
  fgets(line,sizeof line,fin);
  sscanf(line,"%s",filename);
  fgets(line,sizeof line,fin);
  sscanf(line,"%lf %lf %lf %lf %lf %lf",&aM,&eM,&xiM,&OmM,&wM,&mp);
  fgets(line,sizeof line,fin);
  sscanf(line,"%d",&n);
  fgets(line,sizeof line,fin);
  sscanf(line,"%lf %lf %lf %lf %lf %lf",&ip,&a0,&e0,&xi0,&Om0,&w0);
  fclose(fin);

  wM=(wM+OmM)*deg;
  wM=fmod(wM,2.0*pi);
  fprintf(stdout,"wM:%e\n",wM);

  double Tp = pow(a0,1.5);
  double dt = Tp/100.0;
  double t = 0.0;

  fprintf(stdout,"Tp,dt:%e,%e\n",Tp,dt);

  double x[6],xp[6];

  x[1] = a0;
  x[2] = e0;
  x[3] = xi0*deg;
  x[4] = Om0*deg;
  x[5] = (w0+Om0)*deg;
  x[5] = fmod(x[5],2.0*pi);
  fprintf_vec(stdout,"%e ",x,6);

  coefcs(aM,a0);
  
  int jj;
  double te;
  double temp[6],t1;

  gsl_odeiv2_system sys={fcn,NULL,6,NULL};
  gsl_odeiv2_driver* d=
    gsl_odeiv2_driver_alloc_y_new(&sys,
				  gsl_odeiv2_step_rk2,
				  dt,
				  /*epsabs,epsrel*/0.0,1E-6);
  FILE *fl=fopen("results.dat","w");
  int il=0;
  do{
    //fprintf(stdout,"t = %e\n",t);
    gsl_odeiv2_driver_apply_fixed_step(d,&t,dt,1,x);

    for(int k=3;k<=5;k++){
      jj = k-2;
      temp[jj] = x[k]/deg;
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
  }while(t<tfin);
  fclose(fl);
  return 0;
}
