#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutild.h"

//#define NORD 2
#define LEN      1000
//#define EMIN     250.
//#define EMAX     280.
#define NBIN     300

void funcs(double x,double *p,int ma);

int main (void) {

  char input[50],*c,buff[LEN],output[50]="mom_out.dat";
  //  char input[50]="test.dat";
  FILE *fin,*fout,*gp;
  //  double de,ene[BIN_ENE],dbin,ein,amp=1.3;
  //  double *eout,*eo,*pout,dE;
  int i,j;
  int xbin,ybin,mat[NBIN][NBIN],nd,*ia,nord;
  double x,y,xmin,xmax,ymin,ymax,*xx,*yy,*sig,*cf,**covar,chisq;


  printf("order of polynomial : ");
  scanf("%d", &nord);

  //  nord=NORD;

  printf("input file : ");
  scanf("%s", &input);

  printf("input Xmin Xmax : ");
  scanf("%lf %lf", &xmin,&xmax);
  printf("input Ymin Ymax : ");
  scanf("%lf %lf", &ymin,&ymax);


  fin=fopen(input, "r");
  if(fin==NULL) {
    printf("can't open %s",input);
    exit(1);
  }
  xbin=ybin=0;
  while(fgets(buff,LEN,fin) != NULL) {
    c=buff;
    c=strtok(c, ",\n \r");
    //    mat[ybin][xbin]=atoi(c);
    mat[xbin][ybin]=atoi(c);
    c=NULL;
    ybin++;
    if(ybin==NBIN) {
      xbin++;
      ybin=0;
    }
  }
  fclose(fin);

  xx=dvector(1,NBIN*NBIN);
  yy=dvector(1,NBIN*NBIN);
  sig=dvector(1,NBIN*NBIN);
  cf=dvector(1,nord);
  ia=ivector(1,nord);
  covar=dmatrix(1,nord,1,nord);

  nd=1;
  fout=fopen(output, "w");
  for(j=0;j<NBIN;j++) {
    for(i=0;i<NBIN;i++) {
      if(mat[j][i]>0) {
	x=i*(xmax-xmin)/NBIN+xmin;
	y=j*(ymax-ymin)/NBIN+ymin;
	xx[nd]=x;
	yy[nd]=y;
	sig[nd]=1./sqrt(mat[j][i]);
	nd++;
	//	fprintf(fout, "%lf %lf %d\n", x,y,mat[i][j]);
	//	fprintf(fout, "%lf  %lf  %d\n", x,y,mat[j][i]);
	fprintf(fout, "%lf  %lf  %lf\n", x,y,1/sqrt(mat[j][i]));
      }
    }
  }
  fclose(fout);

  for(i=1;i<=nord;i++) ia[i]=i;
  lfit(xx,yy,sig,nd-1,cf,ia,nord,covar,&chisq,funcs);

  /* fitting by using gnuplot */
  gp=popen("gnuplot -geometry 800x800 -persist","w");
//  fprintf(gp, "f(x)=p0");
//  for(i=2;i<=nord;i++) {
//    fprintf(gp, "+p%d*x**%d",i-1,i-1);
//  }
//  fprintf(gp, "\n");
//  fprintf(gp, "fit f(x) \"%s\" using 1:2:3 via p0",output);
//  for(i=2;i<=nord;i++) {
//    fprintf(gp, ",p%d", i-1);
//  }
//  fprintf(gp,"\n");
//  fprintf(gp, "plot \"%s\" using 1:2 w p pt 0, f(x)\n",output);
  fprintf(gp, "f(x)=%lf",cf[1]);
  for(i=2;i<=nord;i++) {
    fprintf(gp, "+%lf*x**%d",cf[i],i-1);
  }
  fprintf(gp, "\n");
  fprintf(gp, "plot \"%s\" using 1:2 w p pt 0, f(x)\n",output);
  fclose(gp);

  printf("Coefficient values\n");
  for(i=1;i<=nord;i++) {
    printf("p[%d]=%lf (+- %lf)\n",i-1,cf[i],sqrt(covar[i][i]));
  }
  printf("paste to encpid.f\n");
  for(i=1;i<=nord;i++) {
    printf("       p%d=%.10lf\n",i-1,cf[i]);
  }



  return 0;
}

void funcs(double x,double *p,int ma) {
  int i;
  p[1]=1.0;
  for(i=2;i<=ma;i++) p[i]=p[i-1]*x;
}

