void Drift2(double el, double inangle, double outangle, double vin[2][3], double vout[2][3], double *pl){
  double t;
  double xe, ye, ze;
  double x, y, z;
  double dy, dz;
  //In
  x = vin[0][0];
  y = vin[0][1];
  z = vin[0][2];
  dy= vin[1][1];
  dz= vin[1][2];
  
  //Out
  t = tan(outangle);
  ye = (dy*(el-x)+y)/(1.-t*dy);
  xe = t*ye;
  ze = z+dz*(xe+el-x);
  vout[0][0] = xe;
  vout[0][1] = ye;
  vout[0][2] = ze;
  vout[1][0] = 1.;
  vout[1][1] = dy;
  vout[1][2] = dz;
  *pl +=sqrt(1.+dz*dz)*sqrt((xe+el-x)*(xe+el-x)+(ye-y)*(ye-y)); 
  return;
}
