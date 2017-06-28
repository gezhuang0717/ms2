// replace invalid value of error to zero
// wroted by S.S. (2016/1/9)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LEN      1000
#define NDAT     2000

int main (int argc, char *argv[]) {

  char input[50],*c,buff[LEN],label[100],title[NDAT][100],output[100];
  FILE *fin,*fout;
  int i,j,ndat,id[NDAT];
  double coef[NDAT][3],dcoef[NDAT][3];

   /* read input file */
 // printf("input file <Enter:2dvec.dat>: ");
 // getline(input,50);
  //  if(strcmp(input,"\0")==0) strcpy(input,"2dvec.dat");
  if(argv[1]==NULL) {
    printf("specify input file\n");
    exit(1);
  }
  fin=fopen(argv[1], "r");
  if(fin==NULL) {
    printf("can't open %s\n",argv[1]);
    exit(1);
  }
  strncpy(output,argv[1],50);

  ndat=0;
  while(fgets(buff,LEN,fin) != NULL) {
    c=buff;
    i=0;
    if(ndat==0) {
      strcpy(label,c);
    }
    while((c=strtok(c, ",'\n\r")) != NULL) {
      if(i==0) id[ndat]=atoi(c);;
      if(i==1) coef[ndat][0]=strtod(c,NULL);
      if(i==2) dcoef[ndat][0]=strtod(c,NULL);
      if(i==3) coef[ndat][1]=strtod(c,NULL);
      if(i==4) dcoef[ndat][1]=strtod(c,NULL);
      if(i==5) coef[ndat][2]=strtod(c,NULL);
      if(i==6) dcoef[ndat][2]=strtod(c,NULL);
      if(i==7) strcpy(title[ndat],c);
      c=NULL;
      i++;
    }
    ndat++;
  }
  fclose(fin);

  /* output file for GNUPLOT */
  fout=fopen(output,"w");
  fprintf(fout, "%s", label);
  for(i=1;i<ndat;i++) {
    fprintf(fout, "%d,", id[i]);
    for(j=0;j<3;j++) {
      if(coef[i][0]==0.) fprintf(fout, ",,");
      else fprintf(fout, "%lf,%lf,", coef[i][j],dcoef[i][j]);
    }
    fprintf(fout, "%s,\n", title[i]);
  }
  fclose(fout);

  printf("output file : %s\n", output);

  return 0;
}

//void get_line(char s[], int lim) {
//  int c,i;
//
//  for(i=0; i<lim-1 && (c=getchar()) != '\n'; ++i)
//    s[i] = c;
//
//  s[i] = '\0';
//}
