/* Estimate inelastic events with dp/p histgram [S.Suzuki] */
/* FORTRAN programs are used to particles sim. and maximum likelihood fitting */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutild.h"

#define LEN      1000  /* Number of character */
#define PMAX57   3.    /* Maximum dp/p of corrected delta57 */
#define PMIN57   -3.   /* Minimum dp/p of corrected delta57 */
#define NCF      2     /* Number of fitting parameter */
#define BINMAX   1000  /* Maximum BIN size */
#define FTOL     1.e-4  /* Tollerance in function, amoeba */
#define NEVENT   1000000   /* Number of events in simulation */

void mlfit(double xx[],double yy[],int ndat,double cf[],double dcf[],int nc);
double ran1(long *iseed);
double chifunc(double *p);
void mkhist(double val,double min,double max,int bin,double count[]);
double dp2e(int nmass,int natom,double dp,double brho);
double e2dp(int nmass,int natom,double ein,double brho);
void fitthick(double thickini,double *thick,double *offset);
//void fitthick(double A,double Z,double slope,int isotp,double ex,double ec,double neve,double *thick);
double func(double);

/* shared variables */
double *d57,*cfit,*csim,*csim2,aa[NCF],sa[NCF];
int limit,bin57,np1;
double *eex;
double slope,exex,ein,emin,emax,sig_exp,ave,sig; /* for using in func */
int mass,charge,istp,neve; /* for using in func */

int main (void) {

  char input[50],file57[50],tmp[50],output[50]="sim/inelsim_",*c,buff[LEN],eps[50]="eps/";
  char test[50]="test.dat";
  FILE *fin,*fout,*gp;
  double ec,ex2,dp57,dp,dp35,ccc[BINMAX];
  double BRHO35,BRHO57,thickt,cf1[NCF],dcf1[NCF],cf2[NCF],dcf2[NCF];
  double inel_max,inel_min,ninel,dinel,nexp,rinel,drinel;
  double ex,**p,*y,*min,*max,csqo=0.,chisq,chisq2,csqtgt;
  double sum1,sum2,cmax=0.,amp,prm[2][2],chit,*ct,weight;
  double *cexp,*csim_max,*csim_min,*eout,*pout,*dpcor,pave;
  double eave,*ecor,de,eoc;
  double cexp_max,cexp_l,cexp_h,fwhm_exp,thick_cor,eoffset;
  double csim_max_val,max_csim,fwhm_sim,sig_sim;
  int jl,jh,simjl,simjh;
  int i,j,k,l,mode,ol,nl,ndim,nfunk,jcmax,offset,A,Z,isotp,nevent,np2;
  long iseed={8234587};

  /* select mode */
  printf("input mode# (0 or 1) : ");
  scanf("%d", &mode);
  if(mode!=0 && mode!=1) {
    printf("please input correct number.\n");
    exit(1);
  }

  /* read parameter file */
  printf("input file : ");
  scanf("%s", &input);
  fin=fopen(input, "r");
  if(fin==NULL) {
    printf("can't open %s\n",input);
    exit(1);
  }
  j=0;
  while(fgets(buff,LEN,fin) != NULL) {
    c=buff;
    c=strtok(c, ", \n\r");
    if(j==0) A=atoi(c);
    if(j==1) Z=atoi(c);
    if(j==2) dp35=strtod(c,NULL);
    if(j==3) BRHO35=strtod(c,NULL);
    if(j==4) BRHO57=strtod(c,NULL);
    if(j==5) strcpy(file57,c);
    if(j==6) limit=atoi(c);
    if(j==7) csqtgt=strtod(c,NULL);
    if(j==8) thickt=strtod(c,NULL);
    //    if(j==9) nevent=atoi(c);
    c=NULL;
    j++;
  }
  fclose(fin);  /* end of parameter reading */
  nevent=NEVENT;  /* number of events in simulation */

  if(mode==1) {
    printf("input slope  : ");
    scanf("%lf", &cf1[1]);
//    cf2[1]=3.3132*cf1[1];     /* empirical value of 23Na */
//    dcf2[1]=0.9936;    /* empirical value of 23Na */
//    cf2[1]=3.6582*cf1[1];     /* empirical value of 23Na */
//    dcf2[1]=0.2448;    /* empirical value of 23Na */
  }

  /* read delta57 data*/
  fin=fopen(file57, "r");
  if(fin==NULL) {
    printf("can't open %s\n",file57);
    exit(1);
  }
  j=0;
  nexp=0.;
  bin57=0;
  while(fgets(buff,LEN,fin) != NULL) {
    c=buff;
    c=strtok(c, ", \n\r");
    ccc[bin57]=strtod(c,NULL);
    nexp+=ccc[bin57];
    c=NULL;
    bin57++;
  }
  fclose(fin);   /* end of data reading */

  /* output file name */
  strcpy(tmp,input);
  strtok(tmp,".");
  strcat(output,tmp);
  strcat(output,".csv");
  strcat(eps,tmp);
  strcat(eps,".eps");

  /* difine vector for amoeba */
  ndim=1;
  min=dvector(1,ndim);
  max=dvector(1,ndim);
  min[1]=1.e-6;               /* minimum value of amplification */
  max[1]=3.e-3;               /* maximum value of amplification */
  p=dmatrix(1,ndim+1,1,ndim);
  y=dvector(1,ndim+1);

  ec=dp2e(A,Z,dp35,BRHO35);     /* beam energy [MeV/u] */
  ex=4.436;  /* the first excitation state of 12C [MeV] */
  ex2=ex;   /* for simulation */

  dp57=(PMAX57-PMIN57)/bin57;  /* dp/p per 1bin */
  d57=dvector(0,bin57);
  cexp=dvector(0,bin57);
  eex=dvector(0,bin57);
  cfit=dvector(0,bin57);
  cexp_max=ccc[bin57/2-1]; /* maximum counts of exp. data */

  for(j=0;j<bin57;j++) {
    d57[j]=j*dp57+PMIN57;
    cexp[j]=ccc[j];
    //    cexp[j]=(ccc[j]>0.?ccc[j]-1.:ccc[j]);
    eex[j]=dp2e(A,Z,d57[j],BRHO57)-dp2e(A,Z,0.,BRHO57);
    if(j<bin57/2&&cexp[j]>cexp_max/2.&&cexp[j-1]<cexp_max/2.)  {
      cexp_l=cexp[j];
      jl=j;
    }
    if(j>bin57/2&&cexp[j]<cexp_max/2.&&cexp[j-1]>cexp_max/2.)  {
      cexp_h=cexp[j];
      jh=j;
    }
  }
  /* calculate standerd deviation for exp. data */
  printf("%d,%d\n", jl,jh);
  fwhm_exp=eex[jh]-eex[jl]; /* full width at half maximum */
  sig_exp=fwhm_exp/2./sqrt(2.*log(2));
  printf("%lf, sig(exp)=%lf\n", fwhm_exp,sig_exp);


  if(mode==0) {
    mlfit(eex,cexp,limit,cf1,dcf1,NCF);  /* MLH fitting by fortran pragram */
    mlfit(d57,cexp,limit,cf2,dcf2,NCF);  /* MLH fitting by fortran pragram */
    printf("Maximum likelihood fitting finished!\n");
  }
  if(mode==1) {  /* fixed slope mode */
    /* coefficient for Eex axis */
    sum1=sum2=0.;
    for(j=0;j<=limit;j++) {
      sum1+=cexp[j];
      sum2+=exp(cf1[1]*eex[j]);
    }
    cf1[0]=sum1/sum2;            /* Amplification */
    dcf1[0]=sqrt(sum1/pow(sum2,2));
    dcf1[1]=0.;
//    /* coefficient for delta_57 axis */
//    sum1=sum2=0.;
//    for(j=0;j<=limit;j++) {
//      sum1+=cexp[j];
//      sum2+=exp(cf2[1]*d57[j]);
//    }
//    cf2[0]=sum1/sum2;            /* Amplification */
//    dcf2[0]=sum1/pow(sum2,2);
  }
  for(i=0;i<NCF;i++) printf("cf1[%d] = %lf(+-%lf)\n", i,cf1[i],dcf1[i]);
  //  for(i=0;i<NCF;i++) printf("cf2[%d] = %lf(+-%lf)\n", i,cf2[i],dcf2[i]);
  for(j=0;j<bin57;j++)  cfit[j]=cf1[0]*exp(cf1[1]*eex[j]); /* function value */
  /* secure memory */
  eout=dvector(0,nevent);
  pout=dvector(0,nevent);
  dpcor=dvector(0,nevent);
  csim=dvector(0,bin57);
  csim2=dvector(0,bin57);
  csim_max=dvector(0,bin57);
  csim_min=dvector(0,bin57);
  isotp=1;  /* set isotp=1 to calculate the offset value for dp/p */
  //  sim_(&A,&Z,&isotp,&cf1[1],&thickt,&ex2,&nevent,&ec,eout); /* simulation by fortran program */
//  thick_cor=1.0*thickt;
//  printf("thick=%lf\n", thick_cor);
//  sim_(&A,&Z,&isotp,&cf1[1],&thick_cor,&ex2,&nevent,&ec,eout); /* simulatie Eout fortran */
//  exit(0);

  emin=dp2e(A,Z,PMIN57,BRHO57)-dp2e(A,Z,0.,BRHO57);
  emax=dp2e(A,Z,PMAX57,BRHO57)-dp2e(A,Z,0.,BRHO57);
  printf("emin=%lf, emax=%lf\n", emin,emax);
  eoc=dp2e(A,Z,0.,BRHO57);
  ecor=dvector(0,nevent);

//  if(isotp==1) { /* calculate delta57 offset */
//    eave=0.;
//    for(j=0;j<nevent;j++) eave+=eout[j];
//    eave=eave/nevent;
//  }


  /* calculate target thckness for fitting to exp. data */
//  for(j=0;j<nevent;j++) ecor[j]=eout[j]-eave; /* delta57 correction */
//  for(j=0;j<bin57;j++) csim[j]=0.;    /* initialize vector */
//  for(i=0;i<nevent;i++) mkhist(ecor[i],emin,emax,bin57,csim);   /* make histgram */
//  csim_max_val=csim[bin57/2-1];
//  for(j=0;j<bin57;j++) {
//    if(j<bin57/2&&csim[j]>csim_max_val/2.&&csim[j-1]<csim_max_val/2.) simjl=j;
//    if(j>bin57/2&&csim[j]<csim_max_val/2.&&csim[j-1]>csim_max_val/2.) simjh=j;
//  }
//  fwhm_sim=eex[simjh]-eex[simjl];
//  sig_sim=fwhm_sim/2./sqrt(2.*log(2));
//  printf("sig(sim)=%lf\n", sig_sim);
//  printf("%lf, %d,%d\n", csim_max_val,simjl,simjh);


/* for func */
  mass=A;
  charge=Z;
  slope=cf1[1];
  istp=1;
  exex=ex2;
  ein=ec;
  neve=nevent;

  printf("calculate thickness for fitting to Exp. data...\n");
  fitthick(thickt,&thick_cor,&eoffset);
  printf("corrected thickness = %lf\n", thick_cor);
  printf("energy offset = %lf\n", eoffset);


  /* test plot */
//  sim_(&A,&Z,&isotp,&cf1[1],&thick_cor,&ex2,&nevent,&ec,eout); /* simulatie Eout fortran */
//  if(isotp==1) { /* calculate delta57 offset */
//    eave=0.;
//    for(j=0;j<nevent;j++) eave+=eout[j];
//    eave=eave/nevent;
//  }
////  fout=fopen(test,"w");
////  for(i=0;i<nevent;i++) fprintf(fout,"%lf\n",eout[i]-eave);
////  fclose(fout);
////  exit(0);
//  for(j=0;j<bin57;j++) csim[j]=0.;    /* initialize vector */
//  for(i=0;i<nevent;i++) mkhist(eout[i]-eave,emin,emax,bin57,csim);   /* make histgram */
//  csim_max_val=csim[bin57/2-1];
//  fout=fopen(test,"w");
//  for(i=0;i<bin57;i++) fprintf(fout,"%lf %lf %lf\n",eex[i],cexp[i],csim[i]*cexp_max/csim_max_val);
//  fclose(fout);
//  gp=popen("gnuplot -persist", "w");
//  fprintf(gp, "plot \"%s\" using 1:2 w l, \"%s\" using 1:3 w l\n", test,test);
//  fclose(gp);
//  exit(0);



  isotp=2;   /* set isotp=2 to simulate the inelastic event */
  //  sim_(&A,&Z,&isotp,&cf1[1],&thickt,&ex2,&nevent,&ec,eout); /* simulatie Eout fortran */
  sim_(&A,&Z,&isotp,&cf1[1],&thick_cor,&ex2,&nevent,&ec,eout); /* simulatie Eout fortran */
  //  for(j=0;j<nevent;j++) ecor[j]=eout[j]-eave; /* delta57 correction */
  for(j=0;j<bin57;j++) csim[j]=0.;    /* initialize vector */
  //  for(i=0;i<nevent;i++) mkhist(ecor[i],emin,emax,bin57,csim);   /* make histgram */
  for(i=0;i<nevent;i++) mkhist(eout[i]-eoffset,emin,emax,bin57,csim);   /* make histgram */

  /* start calculation for maximum estimation */
  printf("\nnow calculating...\n");
  nl=0;
  do{
    for(i=2;i<=ndim+1;i++) /* initialize matrix for amoeba */                   
      for(j=1;j<=ndim;j++){
	if(j==i-1) p[i][j]=min[j]+(max[j]-min[j])*ran1(&iseed);
	else p[i][j]=p[1][j];
      }
    for(i=1;i<=ndim+1;i++) y[i]=chifunc(&(p[i][1])-1);

    amoeba(p,y,ndim,FTOL,chifunc,&nfunk);  /* simplex minimization */
      for(i=1,ol=0;i<=ndim;i++) {  /* check if parametes are within range */
      if(p[1][i]<min[i] || p[1][i]>max[i]) ol=i;
      }
      printf("Loop#%d: csq=%.1e ",++nl,y[1]);
      printf("%s\n",(ol!=0)?" [out of range]":"");
      if(ol!=0||fabs((y[1]-csqo)/y[1])<1.0E-2||y[1]==0.) {/* if over-range or no change in csq */
	for(i=1;i<=ndim;i++) {
	  printf(" p[%d]=%s%f (min:%s%f,max:%s%f)\n", i,(p[1][i]<0)?"":" ",
		 p[1][i],(min[i]<0)?"":" ",min[i],(max[i]<0) ?"":" ",max[i]);
	  p[1][i]=min[i]+(max[i]-min[i])*ran1(&iseed); /* initialize parameters */
	}
      }
      if(ol!=0) y[1]=1000.;   /* over range */
    csqo=y[1];  /* store chi-square value */
  }while(csqo>csqtgt); /* exit loop */

  chisq=chifunc(p[1]);  /* final chi-square */
  printf("Maximum estimation end!\n");
  printf("chisq=%lf [target=%lf]\n", chisq,csqtgt); 
  printf("cf[0]=%lf, cf[1]=%lf\n", aa[0],aa[1]);

  prm[0][0]=p[1][1];  /* amplification of maximum estimation */
  for(j=0;j<bin57;j++) csim_max[j]=csim2[j];  /* counts of maximum estimation */
  inel_max=0.;
  for(j=0;j<bin57;j++) {
    inel_max+=csim_max[j];   /* maximum events */
  }

//  /* test plot */
//  strcpy(output,"test.csv");
//  fout=fopen(output, "w");
//  for(j=0;j<bin57;j++) {
//    fprintf(fout, "%lf,%lf,%lf,%lf\n",eex[j],cexp[j],sqrt(cexp[j]),csim2[j]);
//  }
//  fclose(fout);
//  /* gnuplot */
//  dp=limit*(PMAX57-PMIN57)/bin57+PMIN57;    /* border line */
//  de=dp2e(A,Z,dp,BRHO57)-eoc;    /* border line */
//  gp=popen("gnuplot -persist","w");
//  fprintf(gp, "set datafile separator \",\"\n");
//  fprintf(gp, "set title \"%s\"\n",tmp);
//  fprintf(gp, "set xlabel \"Earb [MeV/u]\"\n");
//  fprintf(gp, "set ylabel \"Counts\"\n");
//  fprintf(gp, "set logscale y\n");
//  fprintf(gp, "set xrange [-10:10]\n");
//  fprintf(gp, "set yrange [0.1:100000]\n");
//  fprintf(gp, "set arrow from %lf,0.1 to %lf,100000 nohead lc 1\n", de,de);
//  fprintf(gp, "plot \"%s\" using 1:2:3 w errorbars lc 0 title \"Exp. data\",", output);
//  fprintf(gp, "\"%s\" using 1:4 w l lt 1 lc 3 title \"Max. estimation\",", output);
//  fprintf(gp, "%lf*exp(%lf*x) w l title \"MLHfitting Exp.\",",cf1[0],cf1[1]);
//  fprintf(gp, "%lf*exp(%lf*x) w l lc 1 title \"MLHfitting Sim\"\n",aa[0],aa[1]);
//
//  fclose(gp);
//  exit(0);


  /* calculate the minimum value */
  for(j=0;j<bin57;j++) {
    if(csim_max[j]>cmax) {
      cmax=csim_max[j];    /* peak of simulation  */
      jcmax=j;
    }
  }
  amp=p[1][1];   /* amplification of maximum estimation */
  p[1][1]=amp*csim_max[limit]/csim_max[jcmax]; /* amplification of minimum estimation */

  ex2=ex+A*fabs(dp2e(A,Z,d57[jcmax],BRHO57)-dp2e(A,Z,d57[limit],BRHO57)); /* Ex at tail top */
  //  sim_(&A,&Z,&isotp,&cf1[1],&thickt,&ex2,&nevent,&ec,eout); /* simulate Eout by fortran */
  sim_(&A,&Z,&isotp,&cf1[1],&thick_cor,&ex2,&nevent,&ec,eout); /* simulate Eout by fortran */
  for(j=0;j<nevent;j++) ecor[j]=eout[j]-eoffset; /* delta57 correction */
  for(j=0;j<bin57;j++) csim[j]=0.;   /* initialize vector */
  for(i=0;i<nevent;i++) mkhist(ecor[i],emin,emax,bin57,csim);   /* make histgram */
  for(j=0;j<bin57;j++) csim2[j]=p[1][1]*csim[j];

  ct=dvector(0,bin57);
  chisq2=1000.;
  k=0;
  while(k<bin57) {
    for(i=0;i+k<bin57;i++) ct[i+k]=csim2[i];  /* offset counts toword plus axis */
    chit=0.;
    l=0;
    for(j=0;j<=limit;j++) {
      if(cfit[j]>0.) {
	weight=1/sqrt(cfit[j]);
	chit+=pow((cfit[j]-ct[j])*weight,2);
	l++;
      }
    }
    chit=chit/l;
    for(i=0;i<bin57;i++) ct[i]=0.;
    if(chit<=chisq2) {
      chisq2=chit;
      offset=k;
      np2=l;
    }
    k++;
  }

  for(j=0;j<bin57;j++) csim_min[j]=0.;   /* initialize vector */
  for(j=0;j<bin57-offset;j++) csim_min[j+offset]=csim2[j];   /* counts of minimum estimation */

  prm[1][0]=p[1][1]; /* amplification of minimum estimation */
  inel_min=0.;
  for(j=0;j<bin57;j++)  inel_min+=csim_min[j]; /* minimum estimation */

  ninel=(inel_max + inel_min)/2.;    /* average of events */
  dinel=(inel_max-inel_min)/sqrt(12.);   /* SD of continuous uniform distribution */
  rinel=ninel/nexp;                  /* ratio of inel. events to total events */
  drinel=dinel/nexp;                 /* error of ratio */

  /* output results to terminal */
  printf("\n------parameter list-----\n");
  printf("A  Z dp_35  Brho_35  Brho_57  Thickness\n");
  printf("%d %d %lf %lf %lf %lf\n", A,Z,dp35,BRHO35,BRHO57,thickt);
  printf("delta57 file : %s [BIN=%d]\n\n", file57,bin57);
  printf("------simulation results------\n");
  printf("chi-square [Max.Est.] = %.2e, points=%d\n",  chisq,np1);
  printf("chi-square [Min.Est.] = %.2e, points=%d\n",  chisq2,np2);
  printf("amp=%lf\n", amp);
  printf("maximum estimation = %lf\n", inel_max);
  printf("minimum estimation = %lf\n", inel_min);
  printf("Ninel = %lf (+- %lf)\n", ninel,dinel);
  printf("Ninel/Nexp = %.8lf (+-%.8lf)\n",rinel,drinel);

  /* output results to file */
  printf("output file name : %s\n", output);
  fout=fopen(output,"w");
  fprintf(fout, "input file,%s\n", input);
  fprintf(fout, ",A,Z,delta35[%],Brho35,Brho57,Exp. data,BIN,BIN limit,thickness[mm],events\n");
  fprintf(fout, ",%d,%d,%lf,%lf,%lf,%s,%d,%d,%lf,%d\n", A,Z,dp35,BRHO35,BRHO57,file57,bin57,limit,thickt,nevent);
  fprintf(fout, "fitting mode,%d\n",mode);
  fprintf(fout, "coefficient,a1,da1,a2,da2\n");
  fprintf(fout, "Eex axis,%lf,%lf,%lf,%lf\n", cf1[0],dcf1[0],cf1[1],dcf1[1]);
  //  fprintf(fout, "dp/p axis,%lf,%lf,%lf,%lf\n", cf2[0],dcf2[0],cf2[1],dcf2[1]);
  fprintf(fout, "Est.Result,chisq,Nfree,target chisq,amplification,aa[0],sa[0],aa[1],sa[1]\n");
  fprintf(fout, "Max.Est.,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf\n",chisq,np1,csqtgt,prm[0][0],aa[0],sa[0],aa[1],sa[1]);
  fprintf(fout, "Min.Est.,%lf,%d,%lf,%lf\n",chisq2,np2,csqtgt,prm[1][0]);
  fprintf(fout, "Nmax,Nmin,Nave,dN,Nexp,Ninel/Nexp,d(Nexp/Ninel)\n");
  fprintf(fout, "%lf,%lf,%lf,%lf,%.0lf,%.8lf,%.8lf\n\n", 
          inel_max,inel_min,ninel,dinel,nexp,rinel,drinel);
  fprintf(fout, "delta57[%],Exp. counts,Error,Max. counts,Min. counts\n");
  for(j=0;j<bin57;j++) {
    fprintf(fout,"%lf,%lf,%lf,%lf,%lf\n", eex[j],cexp[j],sqrt(cexp[j]),csim_max[j],csim_min[j]);
  }
  fclose(fout);

  printf("output file : %s\n", output);

  /* gnuplot */
  dp=limit*(PMAX57-PMIN57)/bin57+PMIN57;  /* border line */
  de=dp2e(A,Z,dp,BRHO57)-eoc;    /* border line */
  gp=popen("gnuplot -persist","w");
  fprintf(gp, "set datafile separator \",\"\n");
  fprintf(gp, "set title \"%s\"\n",tmp);
  fprintf(gp, "set xlabel \"Earb [MeV/u]\"\n");
  fprintf(gp, "set ylabel \"Counts\"\n");
  fprintf(gp, "set logscale y\n");
  fprintf(gp, "set xrange [-10:10]\n");
  fprintf(gp, "set yrange [0.1:100000]\n");
  fprintf(gp, "set arrow from %lf,0.1 to %lf,100000 nohead lc 1\n", de,de);
  fprintf(gp, "plot \"%s\" using 1:2:3 every ::8 w errorbars lc 0 title \"Exp. data\",", output);
  fprintf(gp, "\"%s\" using 1:4 every ::7 w l lt 1 lc 3 title \"Max. estimation\",", output);
  fprintf(gp, "\"%s\" using 1:5 every ::7 w l lt 1 lc 2 title \"Min. estimation\",",output);
  //  fprintf(gp, "%lf*exp(%lf*x) w l title \"MLHfitting\"\n",cf2[0],cf2[1]);
  fprintf(gp, "%lf*exp(%lf*x) w l title \"MLHfitting\"\n",cf1[0],cf1[1]);
  fprintf(gp, "set output \"%s\"\n", eps);
  fprintf(gp, "set terminal postscript eps enhanced color \"Times-Roman,14\"\n");
  fprintf(gp, "replot\n");
  fprintf(gp, "set term X11\n");
  fclose(gp);

  printf("simulation finished!!\n");
  return 0;
}

void mlfit(double xx[],double yy[],int ndat,double cf[],double dcf[],int nc) {

  int i;
  double *cc,*dc;

  cc=dvector(0,nc);
  dc=dvector(0,nc);
  for(i=0;i<nc;i++) cc[i]=1.;
  mlhfit_(xx,yy,&ndat,cc,dc,&nc);  /* Maximum likelihood fitting by fortran pragram */
  for(i=0;i<nc;i++) {
    cf[i]=cc[i];
    dcf[i]=dc[i];
  }
  free_dvector(cc,0,nc);
  free_dvector(dc,0,nc);
  return;
}

double chifunc(double *p) {
/* return chi-square for maximum estimation */
  double chi2,weight;
  int i,j;

  for(j=0;j<bin57;j++) csim2[j]=p[1]*csim[j];
  //  mlfit(d57,csim2,limit,aa,sa,NCF);  /* MLH fitting by fortran pragram */
  mlfit(eex,csim2,limit,aa,sa,NCF);  /* MLH fitting by fortran pragram */
  chi2=0.;
  np1=0;
  for(j=0;j<=limit;j++) {
    if(cfit[j]>0.) {
      weight=1/sqrt(cfit[j]);
      chi2+=pow((cfit[j]-csim2[j])*weight,2);
      np1++;
    }
  }
  chi2=chi2/(np1-2); /* reduced chi-square devided by degree of feedom */

  return chi2;
}

void mkhist(double val,double min,double max,int bin,double count[]) {
/* make histgram */
  double ch;
  int i;

  ch=(val-min)/(max-min)*bin;
  i=(int)ch;
  if(i>=0 && i<=bin) count[i]+=1.;

  return;
}

double dp2e(int nmass,int natom,double dp,double brho) {
/* convert to dp/p[%] to E[MeV/u] */
  double mu=1.66053892e-27; /* kg */
  double c=2.99792458e+8;   /* m/s */
  double e=1.60217657e-19;  /* C */
  double mo,p,etot,po,ke;

  mo=mu*pow(c,2)/e*1.e-6;   /* ~931.5 [MeV/c2] */
  po=natom*brho*c*1.e-6; /* MeV/c */
  p=po*(1.+dp/100.);     /* MeV/c */
  ke=pow(p,2)+pow(nmass*mo,2);
  ke=sqrt(ke);
  ke=ke-nmass*mo;

  return ke/nmass; /* Energy [MeV/u] */
}

double e2dp(int nmass,int natom,double ein,double brho) {
/* return dp/p, ein [MeV/u] */
  double mu=1.66053892e-27; /* kg */
  double c=2.99792458e+8;   /* m/s */
  double e=1.60217657e-19;  /* C */
  double mo,p,etot,po;

  mo=mu*pow(c,2)/e*1.e-6;   /* ~931.5 [MeV/c2] */
  etot=nmass*(ein+mo);   /* MeV */
  p=sqrt(pow(etot,2)-pow(nmass*mo,2));  /* MeV/c */
  po=natom*brho*c*1.e-6; /* MeV/c */

  return (p-po)/po*100.; /* dP/P [%] */
}

void fitthick(double thickini,double *thick,double *offset) {
  double ax,bx,cx,fa,fb,fc,tol,xmin;
  bx=thickini;
  ax=0.5*bx;
  cx=2.*bx;
  tol=1.e-4;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,func);
  golden(ax,bx,cx,func,tol,&xmin);
  *(thick)=xmin;  /* corrected thickness */
  printf("sig(sim)=%lf\n", sig);
  func(xmin);
  *(offset)=ave;  /* offset of Energy */

  return;
}

double func(double xx) {
/* return difference of standerd deviation */
  double *eout,*counts,cmax,thick;
  int i,j,jl,jh;
  thick=fabs(xx);
  eout=dvector(0,neve);
  counts=dvector(0,bin57);
  sim_(&mass,&charge,&istp,&slope,&thick,&exex,&neve,&ein,eout); /* simulation by fortran program */
  ave=0.;
  for(j=0;j<neve;j++) ave+=eout[j];
  ave=ave/neve; /* offset */
  for(j=0;j<bin57;j++) counts[j]=0.;    /* initialize vector */
  for(j=0;j<neve;j++) mkhist(eout[j]-ave,emin,emax,bin57,counts);   /* make histgram */
  cmax=counts[bin57/2-1]; /* maximum counts */

  for(j=0;j<bin57;j++) {
    if(j<bin57/2 && counts[j]>=cmax/2. && counts[j-1]<cmax/2.) jl=j;
    if(j>bin57/2 && counts[j]<=cmax/2. && counts[j-1]>cmax/2.) jh=j;
  }

  sig=(eex[jh]-eex[jl])/2./sqrt(2.*log(2));
  //  printf("sig=%lf,thick=%lf\n",sig,thick);

  return pow(sig-sig_exp,2);
}

