int DipoleE(struct dipole d, double dp, double vin[2][3], double vout[2][3],
            int k0, int k1,
	    double  Bxa[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxbx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxcx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxdx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxby[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxcy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxdy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxbz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxcz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bxdz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double  Bya[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bybx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bycx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bydx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Byby[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bycy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bydy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bybz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bycz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bydz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double  Bza[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzbx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzcx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzdx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzby[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzcy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzdy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzbz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzcz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double Bzdz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
	    double *pl){
  double rad,gam,r0,phi,ko;
  double er0,dphi;
  double x,y,z,dx,dy,dz;
  double xe,ye,ze,xx,yx,zx;
  double xe1,ye1,ze1,dx1,dy1,dz1,ko1;
  double xe2,ye2,ze2,dx2,dy2,dz2,ko2;
  double xe3,ye3,ze3,dx3,dy3,dz3,ko3;
  double xe4,ye4,ze4,dx4,dy4,dz4,ko4;
  double xt,yt,zt;
  double x0,y0,z0;
  double dsn,dcs,dtn;
  double cb,cc,af;
  double p,q,r;
  double Eq_y();
  int i,j,k;
  void GetExitBend();
  
  //Magnet Parameters
  r0  =d.e_rad;
  phi =d.e_phi;
  x0  = 0.;
  y0  =-r0;
  z0  = 0.;
  
  //Orbit Parameters
  rad=AMU*BC*GAM/QMR/Cus/BTOP;
  rad=rad*(1.+dp)/d.bendfactor;
  gam=GAM*(1.+BC*BC*dp);
  
  //Iteration N_SLICE times
  /*Rotation Matrices */ 
  dphi=phi/210.;
  dsn =sin(dphi);
  dcs =cos(dphi);
  dtn =tan(dphi);
  
  //Beam Parameters
  x =vin[0][0];
  y =vin[0][1];
  z =vin[0][2];
  dx=vin[1][0];
  dy=vin[1][1];
  dz=vin[1][2];
  
  /* Exit Reference Point (xx,yx) */
  er0=r0;
  xx = er0*dsn;
  yx =-er0*(1.-dcs);
  zx = 0;
  
//  //Transformation into In-Coordinate System
//  /* Position */
//  x = x;
//  y = y;
//  z = z;
//  /* Direction */
//  dx =1.;
//  dy =dy;
//  dz =dz;
//  /* Rotation Center */
//  x0= x0;
//  y0= y0;
//  z0= z0;
  
  for(k=k0;k<k1;k++){
    
    // 4th Order Runge-Kutta Method
    /* RK 1st step */
    xt=x; yt=y; zt=z; dx1=dx; dy1=dy; dz1=dz;
    GetExitBend(k, x, y, z, xt, yt, zt, r0, rad, dphi, dcs, dsn, dtn,
		Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy, Bxbz, Bxcz, Bxdz,
		Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy, Bybz, Bycz, Bydz,
		Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy, Bzbz, Bzcz, Bzdz,
                &xe1, &ye1, &ze1, &dx1, &dy1, &dz1, &ko1);
    
    /* RK 2nd step */
    xt=0.5*(x+xe1); yt=0.5*(y+ye1); zt=0.5*(z+ze1); dx2=dx; dy2=dy; dz2=dz;
    GetExitBend(k, x, y, z, xt, yt, zt, r0, rad, dphi, dcs, dsn, dtn,
		Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy, Bxbz, Bxcz, Bxdz,
		Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy, Bybz, Bycz, Bydz,
		Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy, Bzbz, Bzcz, Bzdz,
                &xe2, &ye2, &ze2, &dx2, &dy2, &dz2, &ko2);
    
    /* RK 3rd step */
    xt=0.5*(x+xe2); yt=0.5*(y+ye2); zt=0.5*(z+ze2); dx3=dx; dy3=dy; dz3=dz;
    GetExitBend(k, x, y, z, xt, yt, zt, r0, rad, dphi, dcs, dsn, dtn,
		Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy, Bxbz, Bxcz, Bxdz,
		Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy, Bybz, Bycz, Bydz,
		Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy, Bzbz, Bzcz, Bzdz,
                &xe3, &ye3, &ze3, &dx3, &dy3, &dz3, &ko3);
    
    /* RK 4th step */
    xt=xe3; yt=ye3; zt=ze3; dx4=dx; dy4=dy; dz4=dz;
    GetExitBend(k, x, y, z, xt, yt, zt, r0, rad, dphi, dcs, dsn, dtn,
		Bxa, Bxbx, Bxcx, Bxdx, Bxby, Bxcy, Bxdy, Bxbz, Bxcz, Bxdz,
		Bya, Bybx, Bycx, Bydx, Byby, Bycy, Bydy, Bybz, Bycz, Bydz,
		Bza, Bzbx, Bzcx, Bzdx, Bzby, Bzcy, Bzdy, Bzbz, Bzcz, Bzdz,
                &xe4, &ye4, &ze4, &dx4, &dy4, &dz4, &ko4);
    
    //Get Exit Position/Direction(Final)
    xe=(xe1+2.*xe2+2.*xe3+xe4)/6.;
    ye=(ye1+2.*ye2+2.*ye3+ye4)/6.;
    ze=(ze1+2.*ze2+2.*ze3+ze4)/6.;
    dx=(dx1+2.*dx2+2.*dx3+dx4)/6.;
    dy=(dy1+2.*dy2+2.*dy3+dy4)/6.;
    dz=(dz1+2.*dz2+2.*dz3+dz4)/6.;
    ko=(ko1+2.*ko2+2.*ko3+ko4)/6.;
 
    //Transformation into Out-Coordinate System
    /* Exit Position */
    x=xe-xx;
    y=ye-yx;
    z=ze-zx;
    p=dcs*x-dsn*y;
    y=dsn*x+dcs*y;
    x=p;
    z=z;
    /* Exit Direction */
    p  =dcs*dx-dsn*dy;
    dy =dsn*dx+dcs*dy;
    dx =1.;
    dy/=p;
    dz/=p;
    /* Rotation Center */
    x0-=xx;
    y0-=yx;
    z0-=zx;
    p =dcs*x0-dsn*y0;
    y0=dsn*x0+dcs*y0;
    x0=p;
    z0=z0;
    
    /* Accumulate Path */
    *pl+=ko;
  }
  
//  //Transformation into Standard Coordinate System
//  /* Exit Position */
//  x = x;
//  y = y;
//  z = z;
//  /* Exit Direction */
//  dx = 1.;
//  dy = dy;
//  dz = dz;
  
  //Set Output
  vout[0][0]=x;
  vout[0][1]=y;
  vout[0][2]=z;
  vout[1][0]=dx;
  vout[1][1]=dy;
  vout[1][2]=dz;
  return 0;
}

void GetExitBend(int k, double x, double y, double z, 
		 double xt, double yt, double zt, double r0,
		 double rad, double dphi, double dcs, double dsn, double dtn,
		 double  Bxa[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxbx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxcx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxdx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxby[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxcy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxdy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxbz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxcz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bxdz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double  Bya[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bybx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bycx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bydx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Byby[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bycy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bydy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bybz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bycz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bydz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double  Bza[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzbx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzcx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzdx[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzby[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzcy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzdy[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzbz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzcz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double Bzdz[Z_SLICE+1][R_SLICE+1][N_SLICE+1],
		 double *xe, double *ye, double *ze,
		 double *dx, double *dy, double *dz, double *ko){
  int i, j;
  double qx[4], qy[4], qz[4];
  int l;
  double sph[4], spa[4], spl[4], spm[4], spz[4], spb[4], spc[4], spd[4];
  double alfa, px, py, pz, p, p0, q, q0, r, rad_t;
  double fx,fy,fz,xc,yc,zc;
  double gx,gy,gz,radk,ca,cb,cc;

  // Get Curvature at (xt, yt, zt)
  alfa = atan(xt/(yt+r0));
  r = alfa/dphi;
  q = sqrt(xt*xt+(yt+r0)*(yt+r0))-r0;
  j = (int)(q/0.005+16.); q0 = ((double)((int)(q/0.005+16.)))*0.005-0.08;
  p = zt;
  i = (int)(p/0.02+2.); p0 = ((double)((int)(p/0.02+2.)))*0.02-0.04;

  qx[0] = Bxa[i][j-1][k]+Bxbx[i][j-1][k]*r+Bxcx[i][j-1][k]*r*r
         +Bxdx[i][j-1][k]*r*r*r;
  qx[1] = Bxa[i][j][k]  +Bxbx[i][j][k]*r  +Bxcx[i][j][k]*r*r
         +Bxdx[i][j][k]*r*r*r;
  qx[2] = Bxa[i][j+1][k]+Bxbx[i][j+1][k]*r+Bxcx[i][j+1][k]*r*r
         +Bxdx[i][j+1][k]*r*r*r;
  qx[3] = Bxa[i][j+2][k]+Bxbx[i][j+2][k]*r+Bxcx[i][j+2][k]*r*r
         +Bxdx[i][j+2][k]*r*r*r;
  for(l=0;l<3;l++){
      sph[l] = 0.005;
  }
  for(l=1;l<3;l++){
      spa[l] = 3./sph[l]*(qx[l+1]-qx[l])-3./sph[l-1]*(qx[l]-qx[l-1]);
  }
  spl[0] = 1.; spm[0] = 0.; spz[0] = 0.;
  for(l=1;l<3;l++){
      spl[l] = 2.*(sph[l]+sph[l-1])-sph[l-1]*spm[l-1];
      spm[l] = sph[l]/spl[l];
      spz[l] = (spa[l]-sph[l-1]*spz[l-1])/spl[l];
  }
  spl[3] = 1.; spz[3] = 0.; spc[3] = 0;
  for(l=2;l>-1;l--){
      spc[l] = spz[l]-spm[l]*spc[l+1];
      spb[l] = (qx[l+1]-qx[l])/sph[l]-sph[l]*(spc[l+1]+2.*spc[l])/3.;
      spd[l] = (spc[l+1]-spc[l])/3./sph[l];
  }
  px = qx[1]+spb[1]*(q-q0)+spc[1]*(q-q0)*(q-q0)
      +spd[1]*(q-q0)*(q-q0)*(q-q0);
  qy[0] = Bya[i][j-1][k]+Bybx[i][j-1][k]*r+Bycx[i][j-1][k]*r*r
         +Bydx[i][j-1][k]*r*r*r;
  qy[1] = Bya[i][j][k]  +Bybx[i][j][k]*r  +Bycx[i][j][k]*r*r
         +Bydx[i][j][k]*r*r*r;
  qy[2] = Bya[i][j+1][k]+Bybx[i][j+1][k]*r+Bycx[i][j+1][k]*r*r
         +Bydx[i][j+1][k]*r*r*r;
  qy[3] = Bya[i][j+2][k]+Bybx[i][j+2][k]*r+Bycx[i][j+2][k]*r*r
         +Bydx[i][j+2][k]*r*r*r;
  for(l=0;l<3;l++){
      sph[l] = 0.005;
  }
  for(l=1;l<3;l++){
      spa[l] = 3./sph[l]*(qy[l+1]-qy[l])-3./sph[l-1]*(qy[l]-qy[l-1]);
  }
  spl[0] = 1.; spm[0] = 0.; spz[0] = 0.;
  for(l=1;l<3;l++){
      spl[l] = 2.*(sph[l]+sph[l-1])-sph[l-1]*spm[l-1];
      spm[l] = sph[l]/spl[l];
      spz[l] = (spa[l]-sph[l-1]*spz[l-1])/spl[l];
  }
  spl[3] = 1.; spz[3] = 0.; spc[3] = 0;
  for(l=2;l>-1;l--){
      spc[l] = spz[l]-spm[l]*spc[l+1];
      spb[l] = (qy[l+1]-qy[l])/sph[l]-sph[l]*(spc[l+1]+2.*spc[l])/3.;
      spd[l] = (spc[l+1]-spc[l])/3./sph[l];
  }
  py = qy[1]+spb[1]*(q-q0)+spc[1]*(q-q0)*(q-q0)
      +spd[1]*(q-q0)*(q-q0)*(q-q0);
  qz[0] = Bza[i][j-1][k]+Bzbx[i][j-1][k]*r+Bzcx[i][j-1][k]*r*r
         +Bzdx[i][j-1][k]*r*r*r;
  qz[1] = Bza[i][j][k]  +Bzbx[i][j][k]*r  +Bzcx[i][j][k]*r*r
         +Bzdx[i][j][k]*r*r*r;
  qz[2] = Bza[i][j+1][k]+Bzbx[i][j+1][k]*r+Bzcx[i][j+1][k]*r*r
         +Bzdx[i][j+1][k]*r*r*r;
  qz[3] = Bza[i][j+2][k]+Bzbx[i][j+2][k]*r+Bzcx[i][j+2][k]*r*r
         +Bzdx[i][j+2][k]*r*r*r;
  for(l=0;l<3;l++){
      sph[l] = 0.005;
  }
  for(l=1;l<3;l++){
      spa[l] = 3./sph[l]*(qz[l+1]-qz[l])-3./sph[l-1]*(qz[l]-qz[l-1]);
  }
  spl[0] = 1.; spm[0] = 0.; spz[0] = 0.;
  for(l=1;l<3;l++){
      spl[l] = 2.*(sph[l]+sph[l-1])-sph[l-1]*spm[l-1];
      spm[l] = sph[l]/spl[l];
      spz[l] = (spa[l]-sph[l-1]*spz[l-1])/spl[l];
  }
  spl[3] = 1.; spz[3] = 0.; spc[3] = 0;
  for(l=2;l>-1;l--){
      spc[l] = spz[l]-spm[l]*spc[l+1];
      spb[l] = (qz[l+1]-qz[l])/sph[l]-sph[l]*(spc[l+1]+2.*spc[l])/3.;
      spd[l] = (spc[l+1]-spc[l])/3./sph[l];
  }
  pz = qz[1]+spb[1]*(q-q0)+spc[1]*(q-q0)*(q-q0)
      +spd[1]*(q-q0)*(q-q0)*(q-q0);

  fx = 1./sqrt((*dx)*(*dx)+(*dy)*(*dy)+(*dz)*(*dz))*((*dy)*pz-(*dz)*py);
  fy = 1./sqrt((*dx)*(*dx)+(*dy)*(*dy)+(*dz)*(*dz))*((*dz)*px-(*dx)*pz);
  fz = 1./sqrt((*dx)*(*dx)+(*dy)*(*dy)+(*dz)*(*dz))*((*dx)*py-(*dy)*px);
  rad_t = rad/sqrt(fx*fx+fy*fy+fz*fz);

  // Get Center of Motion (xc,yc,zc)
  xc = x+rad_t*fx/sqrt(fx*fx+fy*fy+fz*fz);
  yc = y+rad_t*fy/sqrt(fx*fx+fy*fy+fz*fz);
  zc = z+rad_t*fz/sqrt(fx*fx+fy*fy+fz*fz);

  // Get Exit Position (*xe,*ye,*ze)
  gx = fy*(*dz)-fz*(*dy);
  gy = fz*(*dx)-fx*(*dz);
  gz = fx*(*dy)-fy*(*dx);
  radk = r0+yc-xc/dtn;
  ca = (gx+gy/dtn)*(gx+gy/dtn)+gz*gz/dsn/dsn;
  cb = -2.*radk*(gx*gy+(gy*gy+gz*gz)/dtn);
  cc = (gy*gy+gz*gz)*radk*radk-gz*gz*rad_t*rad_t;
  *xe = (-cb+sqrt(cb*cb-4.*ca*cc))/2./ca+xc;
  *ye = *xe/dtn-r0;
  *ze = -((*xe-xc)*gx+(*ye-yc)*gy)/gz+zc;
  *ko = 2.*asin(sqrt((*xe-x)*(*xe-x)+(*ye-y)*(*ye-y)
                    +(*ze-z)*(*ze-z))/2./rad_t)*rad_t;

  // Get Exit Direction (dx,dy,dz)
  r = gy*(*ze-zc)-gz*(*ye-yc);
  *dy = gz*(*xe-xc)-gx*(*ze-zc);
  *dz = gx*(*ye-yc)-gy*(*xe-xc);
  *dx = 1;
  *dy/= r;
  *dz/= r;

  return;
}
