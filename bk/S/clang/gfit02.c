// by ss (2015/11/26)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutild.h"

#define LEN      1000
#define NBIN     2000
#define MAXITER 1000
#define CHISTABLE 4
#define MINCHISQ 100

void fgauss(double x, double a[], double *y, double dyda[], int na);

int main (void) {

  char input[50],*c,buff[LEN],output[50]="vfit.dat";
  FILE *fin,*fout,*gp;
  int i,j;
  int bin,nd,*ia,nord=3,iter,itst;
  double vec[NBIN],x,y,xmin,xmax,ymin,ymax,*xx,*yy,*sig,*cf,**covar,**alpha,chisq,ochisq,alamda;


//  printf("order of polynomial : ");
//  scanf("%d", &nord);

  printf("input file : ");
  scanf("%s", &input);

  printf("input Xmin Xmax : ");
  scanf("%lf %lf", &xmin,&xmax);


  fin=fopen(input, "r");
  if(fin==NULL) {
    printf("can't open %s",input);
    exit(1);
  }
  bin=0;
  while(fgets(buff,LEN,fin) != NULL) {
    c=buff;
    c=strtok(c, ",\n \r");
    vec[bin]=strtod(c,NULL);
    c=NULL;
    bin++;
  }
  fclose(fin);


  xx=dvector(1,bin);
  yy=dvector(1,bin);
  yy=&vec[-1];
  sig=dvector(1,bin);
  cf=dvector(1,nord);
  ia=ivector(1,nord);
  covar=dmatrix(1,nord,1,nord);
  alpha=dmatrix(1,nord,1,nord);

  for(i=1;i<=bin;i++) {
    xx[i]=(i-1)*(xmax-xmin)/bin+xmin;
    sig[i]=1./sqrt(yy[i]);
  }
  for(i=1;i<=nord;i++) {
    cf[i]=1.;
    ia[i]=1;
  }
  alamda=-1;
  itst=0;
  for(iter=1;iter<=MAXITER;iter++) {
    printf("iter=%d, chisq=%lf, alamda=%lf\n", iter,chisq,alamda);
    ochisq=chisq;
    mrqmin(xx,yy,sig,bin,cf,ia,nord,covar,alpha,&chisq,fgauss,&alamda);
    if(alamda==0.) break;
    if(chisq>ochisq) itst=0;
    else if(fabs(ochisq-chisq)<0.1) itst++;
    if(itst<CHISTABLE) continue;
    if(ochisq>MINCHISQ) continue;

    alamda=0.0;
    break;
  }
  printf("chisq=%lf\nalamda=%lf\n", chisq,alamda);
  for(i=1;i<=nord;i++) printf("cf[%d] = %lf\n", i,cf[i]);




  fout=fopen(output, "w");
  for(j=1;j<=bin;j++) {
    fprintf(fout, "%lf %lf %lf\n", xx[j],yy[j],sig[j]);
  }
  fclose(fout);

  gp=popen("gnuplot -geometry 800x800 -persist","w");
  fprintf(gp, "f(x)=%lf*exp(-((x-%lf)/%lf)**2.)",cf[1],cf[2],cf[3]);
//  fprintf(gp, "f(x)=%lf",cf[1]);
//  for(i=2;i<=nord;i++) {
//    fprintf(gp, "+%lf*x**%d",cf[i],i-1);
//  }
  fprintf(gp, "\n");
  fprintf(gp, "plot \"%s\" using 1:2 w p pt 0, f(x)\n",output);
  fclose(gp);

  printf("Coefficient values\n");
  for(i=1;i<=nord;i++) {
    printf("p[%d]=%lf (+- %lf)\n",i-1,cf[i],sqrt(covar[i][i]));
  }
  printf("paste to encpid.f\n");
  for(i=1;i<=nord;i++) {
    //   printf("       p%d=%.10lf\n",i-1,cf[i]);
    printf("%.6e\n",cf[i]);
  }
//  printf("paste in anapaw\nfunc/plot ");
//  for(i=1;i<=nord;i++) {
//    printf("%+.5e*x**%d", cf[i],i-1);
//  }
//  printf(" %lf %lf\n",xmin,xmax);


  return 0;
}

//  void funcs(double x, double a[], double *y, double dyda[], int na)
//    {
//      int i;
//      double fac,ex,arg;
//
//      *y=0.0;
//      for (i=1;i<=na-1;i+=3) {
//	arg=(x-a[i+1])/a[i+2];
//	ex=exp(-arg*arg);
//	fac=a[i]*ex*2.0*arg;
//	*y += a[i]*ex;
//	dyda[i]=ex;
//	dyda[i+1]=fac/a[i+2];
//	dyda[i+2]=fac*arg/a[i+2];
//      }
//    }


//void funcs(double x,double *p,int ma) {
//  int i;
//  p[1]=1.0;
//  for(i=2;i<=ma;i++) p[i]=p[i-1]*x;
//}

