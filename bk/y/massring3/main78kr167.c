#include <stdio.h>
#include <math.h>

#include "Ring.c"
#include "Off_y.c"

int main(void){
  double in[15], out[15];
  double param[20];
  char   option[128];
  double dy,dyp;
  double dtime0,dtof;
  double y_max,yp_max;
  double dp,am,bc,gm,pe,ee;
  int i,j,k,m,n,l;
  
  double Off_y();
  int Ring();
  
  // Parameters
  param[0]=6.;      /* Number of Sectors                   */
  param[1]=4.045;   /* Central Orbit: Radius[m]            */
  param[2]=60.3507; /* Central Orbit: Circumference[m]     */
  param[3]=0.5310246;
                    /* Central Orbit:  Beta for 167.8MeV/u   */
  param[4]=36./77.9204;
                    /* Central Orbit: Charge to Mass Ratio */
  param[5]=2.;   /* Number of Turns                     */
  param[6]=1400.;   /* Mag. Sector: No. of Slices (t)      */
  param[12]=32.;    /* Mag. Sector: No. of Slices (r)      */
  param[13]=4.;     /* Mag. Sector: No. of Slices (z)      */
  param[7]=0.04;    /* Mag. Sector: Vertical Spacing[m]    */
  param[8]=0.20;    /* Mag. Sector: Horisontal  Spacing[m] */
  param[14]=1.040539775;/* Mag. Sector: Field at Flat Top      */
  param[9]=1.57751; /* Multipl. Factor for w2 NOT USED     */
  param[10]=1.51769;/* Multipl. Factor for w1 NOT USED     */
  param[11]=1.01;/* Multipl. Factor for Bending radius  */
  for(i=15;i<19;i++)
    param[i]=0.;    /* Reserved                            */
  param[19]=1.;     /* Debug: Parameters[2.] Print[1.]     */
		    /*        No Print  [0.]               */
  
  // Standard Tof in use
  dtime0 =param[5]*param[2]/param[3]/2.99792458e8;
  dtime0*=(1.);  /* Fine Adjust */
  
  // Emittance Range
  y_max =0.032;
  yp_max=0.005;
  
  // Momentum
  for(l=0; l<1; l++){
    dp=l*0.001;
    am=AMU;
    bc=param[3];
    gm=1./sqrt(1.-bc*bc);
    pe=am*gm*bc*(1.+dp);
    ee=sqrt(pe*pe+am*am)-am;
    for (m=25; m<26; m++){
      dy = -y_max+2.*y_max*((double)m/50.)+Off_y(dp,param);
      for (n=25; n<26; n++){
	dyp = -yp_max+2.*yp_max*((double)n/50.);
	in[0]=dy;   //y
	in[1]=dyp;  //y'
	in[2]=0.;   //z
	in[3]=0.;   //z'
	in[4]=ee;   //e
	in[5]=0.;   //tof
	in[6]=77.9204;  //m
	in[7]=36.;  //q
	for(i=8;i<15;i++)
	  in[i]=0.;
	for(i=0;i<15;i++)
	  out[i]=0.;
	if(Ring(in, out, param, option)!=-1){
	  dtof=out[11]/dtime0-1.;
	  if(param[19]==0.){
	    printf(" %12.9lf %12.9lf %12.9lf %12e\n", dp, dy, dyp, dtof);
	  }
	}
      }
    }
  }
  return 0;
}
