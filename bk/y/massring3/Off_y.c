double Off_y(double dp, double param[]){
  double phi2;
  double el, del;
  // Ring Parameters
  double  RD, CL, BA, BC;
  int NU;
//  // Set Parameters
//  NU =(int)(param[0]+0.1);     /* Number of Sectors                   */
//  RD =param[1];                /* Central Orbit: Radius[m]            */
//  CL =param[2];                /* Central Orbit: Circumference[m]     */
//  BC =param[3];                /* Central Orbit: Beta                 */
//  BA=2.*Pi/(double)NU;         /* Mag. Sector: Bending Angle[rad]     */
//  //Estimate Offset_y
//  phi2=BA/2.;
//  el  =CL/(double)NU-RD*BA;
//  del =el-(RD*BA+el)*BC*BC;
//  del/=(2.*sin(phi2));
  return 6.295*dp;
}
