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

double ran1(long *iseed);
void fgauss(double x, double a[], double *y, double dyda[], int na);

int main (void) {

  char input[50],*c,buff[LEN],output[50]="vfit.dat";
  FILE *fin,*fout,*gp;
  //  long iseed={8234587};
  int i,j,bin,nd,*ia,nord=3,iter,itst,ci,li;
  double vec[NBIN],x,y,xmin,xmax,xw,ymax,fwhm,yamax;
  double *xx,*yy,*sig,*cf,*cf0,**covar,**alpha,chisq,ochisq,alamda;


  printf("input file : ");
  scanf("%s", &input);
  printf("input Xmin Xmax : ");
  scanf("%lf %lf", &xmin,&xmax);

  /* open input file */
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
  xw=(xmax-xmin)/bin; /* width of x */

  xx=dvector(1,bin);
  yy=dvector(1,bin);
  yy=&vec[-1];
  sig=dvector(1,bin);
  cf=dvector(1,nord);
  cf0=dvector(1,nord);
  ia=ivector(1,nord);
  covar=dmatrix(1,nord,1,nord);
  alpha=dmatrix(1,nord,1,nord);

  for(i=1;i<=bin;i++) {
    xx[i]=(i-1)*(xmax-xmin)/bin+xmin+xw/2.;
    //    sig[i]=1./sqrt(yy[i]);
//    sig[i]=sqrt(yy[i]);
    if(yy[i]>0.) sig[i]=sqrt(yy[i]);
    else sig[i]=1.;
  }
  for(i=1;i<=nord;i++) ia[i]=1;

  /* initial parameter search */
  ymax=0.;
  for(j=1;j<=bin;j++) {
    if(yy[j]>ymax) {
      ymax=yy[j];
      ci=j;
    }
  }
  for(j=1;j<=bin;j++)  if(yy[j]>=ymax/2.) {
    li=j;
    break;
  }
  for(i=0;i<6;i++) {
    if(ymax>pow(10,i)) yamax=pow(10,i+1);
  }
  fwhm=2.*fabs(xx[ci]-xx[li]);
  cf0[1]=yy[ci];
  cf0[2]=xx[ci];
  cf0[3]=fwhm/2./sqrt(log(2)); /* cf[3]=sqrt(2)*sigma */
  for(i=1;i<=nord;i++) cf[i]=cf0[i];

  /* start iteration */
  alamda=-1;
  itst=0;
  mrqmin(xx,yy,sig,bin,cf,ia,nord,covar,alpha,&chisq,fgauss,&alamda);
  for(iter=1;iter<=MAXITER;iter++) {
    ochisq=chisq;
    mrqmin(xx,yy,sig,bin,cf,ia,nord,covar,alpha,&chisq,fgauss,&alamda);
    if(alamda==0.) break;
    if(chisq>ochisq) itst=0;
    else if(fabs(ochisq-chisq)<0.1) itst++;
    if(itst<CHISTABLE) continue;
//    if(fabs(ochisq-chisq)<1.e-3) {
//      for(i=1;i<=nord;i++) cf[i]=cf0[i]*ran1(&iseed);
//      continue;
//    }
    /* finish fitting */
    alamda=0.0;
    mrqmin(xx,yy,sig,bin,cf,ia,nord,covar,alpha,&chisq,fgauss,&alamda);
    break;
  }

  /* fitting result */
  printf("chisq=%lf\nalamda=%lf\n", chisq,alamda);
  //  for(i=1;i<=nord;i++) printf("cf[%d] = %lf\n", i,cf[i]);
  printf("Constant : %lf (+- %.3e)\n", cf[1],sqrt(covar[1][1]));
  printf("Mean : %lf (+- %.3e)\n", cf[2],sqrt(covar[2][2]));
  printf("Sigma : %lf (+- %.3e)\n", cf[3]/sqrt(2),sqrt(covar[3][3]/2.));

  fout=fopen(output, "w");
  for(j=1;j<=bin;j++) {
    fprintf(fout, "%lf %lf %lf\n", xx[j],yy[j],sig[j]);
  }
  fclose(fout);
  gp=popen("gnuplot -geometry 800x800 -persist","w");
  fprintf(gp, "set logscale y\n");
  fprintf(gp, "set yrange [0.1:%lf]\n",yamax);
  fprintf(gp, "set xrange [%lf:%lf]\n",cf[2]-1.,cf[2]+1.);
  fprintf(gp, "f(x)=%lf*exp(-((x-%lf)/%lf)**2)",cf[1],cf[2],cf[3]);
  fprintf(gp, "\n");
  fprintf(gp, "plot \"%s\" using 1:2 w p pt 0, f(x)\n",output);
  fclose(gp);

//  printf("Coefficient values\n");
//  for(i=1;i<=nord;i++) {
//    printf("p[%d]=%lf (+- %.3e)\n",i-1,cf[i],sqrt(covar[i][i]));
//  }


  return 0;
}
