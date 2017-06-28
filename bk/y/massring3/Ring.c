#include <stdio.h>
#include <math.h>

// Constants
#define Pi  3.1415926535   /* Pi                           */
#define Cus 2.99792458e8   /* Light Velocity[m/sec]        */
#define AMU 931.49404283e6 /* Atomic Mass Unit[MeV/c2]     */

// Ring Parameters
double  RD, BA, GAM, QMR, CL, BC, CLT, BTOP;
double  V_space, H_space;
double  Print;
int NU,N;
int N_SLICE, R_SLICE, Z_SLICE;

// Sector Parameters
struct dipole {
  double e_rad;            /* Central Orbit Radius         */
  double e_phi;            /* Bending Angle                */
  double bendfactor;       /* Mag. Field Correction        */
};

int Ring(double *in, double *out, double *param, char *option){
  FILE *fp;
  double vin[2][3], vout[2][3],pl,pl0;
  double e0,e,dpp,ee,bt;
  double mass,charge,scale;
  struct dipole d;
  double phi2;
  double st,st2,delta;
  double dy,dyp,xb;
  double tof,s;
  int i,j,k,m,n;
  double r,q,p;
  char *name[6];
  
  double Eq_y();
  int Dipole6();
  int DipoleE();
  void Drift2();
  void GetSpline();
  
  // Parameters
  NU =(int)(param[0]+0.1);     /* Number of Sectors                   */
  RD =param[1];                /* Central Orbit: Radius[m]            */
  CL =param[2];                /* Central Orbit: Circumference[m]     */
  BC =param[3];                /* Central Orbit: Beta                 */
  GAM=1./sqrt(1.-BC*BC);       /* Central Orbit: Gamma Factor         */
  QMR=param[4];                /* Central Orbit: Charge to Mass Ratio */
  N =(int)(param[5]+0.1);      /* Number of Turns                     */
  CLT=CL*(double)N;            /* Central Orbit: Total Path Length    */
  N_SLICE=(int)(param[6]+0.1); /* Mag. Sector: No. of Slices (t)      */
  R_SLICE=(int)(param[12]+0.1);/* Mag. Sector: No. of Slices (r)      */
  Z_SLICE=(int)(param[13]+0.1);/* Mag. Sector: No. of Slices (z)      */
  V_space=param[7];            /* Mag. Sector: Vertical Spacing[m]    */
  H_space=param[8];            /* Mag. Sector: Horisontal  Spacing[m] */
  BA=2.*Pi/(double)NU/4.;      /* Mag. Sector: Bending Angle[rad]     */
  BTOP=param[14];              /* Mag. Sector: Field at Flat Top      */
  Print  =param[19];           /* Track: No Print Out[0]/Print Out[1] */
  
  k=N_SLICE+1;
  j=R_SLICE+1;
  i=Z_SLICE+1;
  double Bxfield[i][j][k], Byfield[i][j][k], Bzfield[i][j][k];
  double Bxa[i][j][k], Bxbx[i][j][k], Bxcx[i][j][k], Bxdx[i][j][k];
  double               Bxby[i][j][k], Bxcy[i][j][k], Bxdy[i][j][k];
  double               Bxbz[i][j][k], Bxcz[i][j][k], Bxdz[i][j][k];
  double Bya[i][j][k], Bybx[i][j][k], Bycx[i][j][k], Bydx[i][j][k];
  double               Byby[i][j][k], Bycy[i][j][k], Bydy[i][j][k];
  double               Bybz[i][j][k], Bycz[i][j][k], Bydz[i][j][k];
  double Bza[i][j][k], Bzbx[i][j][k], Bzcx[i][j][k], Bzdx[i][j][k];
  double               Bzby[i][j][k], Bzcy[i][j][k], Bzdy[i][j][k];
  double               Bzbz[i][j][k], Bzcz[i][j][k], Bzdz[i][j][k];
  name[0] = "field.dat";

  // Default Setup for *out
  for(i=0;i<14;i++)
    out[i] = in[i];
  
  // Beam Parameters
  e0 =AMU*(GAM-1.);
  e  =in[4];
  mass=in[6];
  charge=in[7];
  scale=QMR/(charge/mass);
  dpp=scale*sqrt((2.*AMU+e)*e/(2.*AMU+e0)/e0)-1.;
  
  // Initialize Magnetic Sector
  d.e_rad=RD;
  d.e_phi=BA;
  phi2   =BA/2.;
  st     =CL/(double)NU-RD*BA*4.-0.6-0.6-0.6;
  st2    =st/2.-0.5;
  d.bendfactor = param[11];

  // Print Log
  if(Print==2.){
    printf("Central Orbit Radius:%lf [m]\n",RD);
    printf("Bending Angle/1 Unit:%lf [rad]\n",BA);
    printf("Straight Section    :%lf [m]\n",st);
    printf("Momentum Fraction   :%lf \n", 1.+dpp);
  }
  dy = in[0];
  xb = 0;
  vin[0][0] = xb    ;//x
  vin[1][0] = 1.    ;//x'
  vin[0][1] = in[0] ;//y
  vin[1][1] = in[1] ;//y'
  vin[0][2] = in[2] ;//z
  vin[1][2] = in[3] ;//z'
  pl = 0.;
  s=0.;

  // Scan the Field File
  /* Open the file */
  fp = fopen(name[0], "r");
  /* Scan the file */
  for(i=0;i<Z_SLICE+1;i++){
      for(j=0;j<R_SLICE+1;j++){
	  for(k=0;k<N_SLICE+1;k++){
	      fscanf(fp, "%lf %lf %lf %lf %lf %lf\n",
		     &p, &q, &r,
		     &Bxfield[i][j][k],
		     &Byfield[i][j][k],
		     &Bzfield[i][j][k]);
	  }
      }
  }
  /* Close the file */
  fclose(fp);

  /* Make the continuous field with spline */
  GetSpline(N_SLICE, R_SLICE, Z_SLICE, Bxfield,
	    Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy, Bxbz, Bxcz, Bxdz);
  GetSpline(N_SLICE, R_SLICE, Z_SLICE, Byfield,
	    Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy, Bybz, Bycz, Bydz);
  GetSpline(N_SLICE, R_SLICE, Z_SLICE, Bzfield,
	    Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy, Bzbz, Bzcz, Bzdz);

  // Turn the Ring //
  for(i=0;i<N;i++){
    for(j=0;j<NU;j++){
      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      Drift2(st2-1.35, 0., 0., vin, vout,&pl);
      vin[0][0] = vout[0][0];//x Drift 1.511225 m
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=st2-1.35+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(Dipole6(d, dpp, vin, vout, 0, 100,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x Drift 0.5 m
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=0.5+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(DipoleE(d, dpp, vin, vout, 100, 310,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x Enterance Mag
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=RD*BA+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(Dipole6(d, dpp, vin, vout, 310, 430,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x Drift 0.6 m
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=0.6+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(DipoleE(d, dpp, vin, vout, 430, 640,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x TARN II 1
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=RD*BA+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(Dipole6(d, dpp, vin, vout, 640, 760,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x Drift 0.6 m
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=0.6+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(DipoleE(d, dpp, vin, vout, 760, 970,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x TARN II 2
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=RD*BA+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(Dipole6(d, dpp, vin, vout, 970, 1090,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x Drift 0.6 m
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=0.6+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(DipoleE(d, dpp, vin, vout, 1090, 1300,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x Exit Mag
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=RD*BA+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      if(Dipole6(d, dpp, vin, vout, 1300, 1400,
		 Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy,
                                        Bxbz, Bxcz, Bxdz,
		 Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy,
                                        Bybz, Bycz, Bydz,
		 Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy,
                                        Bzbz, Bzcz, Bzdz,		 
		 &pl)==-1)
	return -1;
      vin[0][0] = vout[0][0];//x Drift 0.5 m
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=0.5+vin[0][0];

      if(Print==1.&&i<10000)
	printf("%14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",s,
	       vin[0][1],vin[1][1],vin[0][2],vin[1][2]);
      Drift2(st2+1.35, 0., 0., vin, vout,&pl);
      vin[0][0] = vout[0][0];//x Drift 1.511225 m
      vin[1][0] = vout[1][0];//x'
      vin[0][1] = vout[0][1];//y
      vin[1][1] = vout[1][1];//y'
      vin[0][2] = vout[0][2];//z
      vin[1][2] = vout[1][2];//z'
      s+=st2+1.35+vin[0][0];
    }
  }
  ee=e+AMU;
  bt=sqrt(ee*ee-AMU*AMU)/ee;

  /* Path Correction */
  pl -=vin[0][0];
  
  tof = pl/bt/Cus;
  out[0] = vin[0][1] ;//y
  out[1] = vin[1][1] ;//y'
  out[2] = vin[0][2] ;//z
  out[3] = vin[1][2] ;//z'
  out[5]+= tof       ;//time
  out[11]= tof       ;//tof
  return 0;
}

#include "DM_slice6.c"
#include "DM_sliceE.c"
#include "ST2.c"
#include "Eq_y.c"

void GetSpline(int N_SLICE, int R_SLICE, int Z_SLICE,
	       double field[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double  a[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double bx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double cx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double ex[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double by[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double cy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double ey[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double bz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double cz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	       double ez[Z_SLICE+1][R_SLICE+1][N_SLICE+1]){
    int i, j, k;
    k=N_SLICE+1;
    double sph[k], spa[k], spl[k], spm[k], spz[k];
    
    // Make Continuous Data
    /* theta(k)-direction r(j)-fix z(i)-fix*/
    for(i=0;i<Z_SLICE+1;i++){
	for(j=0;j<R_SLICE+1;j++){
	    for(k=0;k<N_SLICE+1;k++){
		a[i][j][k] = field[i][j][k];
	    }
	    for(k=0;k<N_SLICE;k++){
		sph[k] = 1.;
	    }
	    for(k=1;k<N_SLICE;k++){
		spa[k] = 3./sph[k]*(a[i][j][k+1]-a[i][j][k])
		        -3./sph[k-1]*(a[i][j][k]-a[i][j][k-1]);
	    }
	    spl[0] = 1.;
	    spm[0] = 0.;
	    spz[0] = 0.;
	    for(k=1;k<N_SLICE;k++){
		spl[k] = 2.*(sph[k]+sph[k-1])-sph[k-1]*spm[k-1];
		spm[k] = sph[k]/spl[k];
		spz[k] = (spa[k]-sph[k-1]*spz[k-1])/spl[k];
	    }
	    spl[N_SLICE] = 1.;
	    spz[N_SLICE] = 0.;
	    cx[i][j][N_SLICE] = 0.;
	    for(k=N_SLICE-1;k>-1;k--){
		cx[i][j][k] = spz[k]-spm[k]*cx[i][j][k+1];
		bx[i][j][k] = (a[i][j][k+1]-a[i][j][k])/sph[k]
		             -sph[k]*(cx[i][j][k+1]+2.*cx[i][j][k])/3.;
		ex[i][j][k] = (cx[i][j][k+1]-cx[i][j][k])/3./sph[k];
	    }
	}
    }

    /* r(j)-direction theta(k)-fix z(i)-fix*/
    for(i=0;i<Z_SLICE+1;i++){
	for(k=0;k<N_SLICE+1;k++){
	    for(j=0;j<R_SLICE;j++){
		sph[j] = 0.005;
	    }
	    for(j=1;j<R_SLICE;j++){
		spa[j] = 3./sph[j]*(a[i][j+1][k]-a[i][j][k])
		        -3./sph[j-1]*(a[i][j][k]-a[i][j-1][k]);
	    }
	    spl[0] = 1.;
	    spm[0] = 0.;
	    spz[0] = 0.;
	    for(j=1;j<R_SLICE;j++){
		spl[j] = 2.*(sph[j]+sph[j-1])-sph[j-1]*spm[j-1];
		spm[j] = sph[j]/spl[j];
		spz[j] = (spa[j]-sph[j-1]*spz[j-1])/spl[j];
	    }
	    spl[R_SLICE] = 1.;
	    spz[R_SLICE] = 0.;
	    cy[i][R_SLICE][k] = 0.;
	    for(j=R_SLICE-1;j>-1;j--){
		cy[i][j][k] = spz[j]-spm[j]*cy[i][j+1][k];
		by[i][j][k] = (a[i][j+1][k]-a[i][j][k])/sph[j]
		             -sph[j]*(cy[i][j+1][k]+2.*cy[i][j][k])/3.;
		ey[i][j][k] = (cy[i][j+1][k]-cy[i][j][k])/3./sph[j];
	    }
	}
    }

    /* z(i)-direction theta(k)-fix r(j)-fix*/
    for(j=0;j<R_SLICE+1;j++){
	for(k=0;k<N_SLICE+1;k++){
	    for(i=0;i<Z_SLICE;i++){
		sph[i] = 0.02;
	    }
	    for(i=1;i<Z_SLICE;i++){
		spa[i] = 3./sph[i]*(a[i+1][j][k]-a[i][j][k])
		        -3./sph[i-1]*(a[i][j][k]-a[i-1][j][k]);
	    }
	    spl[0] = 1.;
	    spm[0] = 0.;
	    spz[0] = 0.;
	    for(i=1;i<Z_SLICE;i++){
		spl[i] = 2.*(sph[i]+sph[i-1])-sph[i-1]*spm[i-1];
		spm[i] = sph[i]/spl[i];
		spz[i] = (spa[i]-sph[i-1]*spz[i-1])/spl[i];
	    }
	    spl[Z_SLICE] = 1.;
	    spz[Z_SLICE] = 0.;
	    cz[Z_SLICE][j][k] = 0.;
	    for(i=Z_SLICE-1;i>-1;i--){
		cz[i][j][k] = spz[i]-spm[i]*cz[i+1][j][k];
		bz[i][j][k] = (a[i+1][j][k]-a[i][j][k])/sph[i]
		             -sph[i]*(cz[i+1][j][k]+2.*cz[i][j][k])/3.;
		ez[i][j][k] = (cz[i+1][j][k]-cz[i][j][k])/3./sph[i];
	    }
	}
    }
    
    return;
}
