#include "fargo3d.h"
			       
double dotproduct(vector v1, vector v2)
{
  double prod;

  prod = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;

  return prod;
}

point pointat(double x, double y, double z)
{
  point pt;
  pt.x = x;
  pt.y = y;
  pt.z = z;
  return pt;
}

point cylin2cart(double x, double phi, double z, double a_au)
{
  point pt;
  double r;
  r = x+a_au;

  pt.x = r*cos(phi);
  pt.y = r*sin(phi);
  pt.z = z;

  return pt;
}

vector vectorbtwnpts(point p1, point p2)
{
  vector v;

  v.x = p2.x - p1.x;
  v.y = p2.y - p1.y;
  v.z = p2.z - p1.z;

  return v;
}

vector crossproduct(vector v1, vector v2)
{
  vector prod;
  /* |  i    j    k   |
   * | v1.x v1.y v1.z |
   * | v2.x v2.y v2.z |
   */ 
  prod.x = v1.y*v2.z - v1.z*v2.y;
  prod.y = v1.z*v2.x - v1.x*v2.z;
  prod.z = v1.x*v2.y - v1.y*v2.x;

  return prod;
}

vector unitvector(vector v)
{
  vector vans;
  double scale;
  scale = sqrt(dotproduct(v,v));
  vans.x = v.x/scale;
  vans.y = v.y/scale;
  vans.z = v.z/scale;
  return vans;
}

double F_incident(double T_star, double r_star, double a_AU)
{
  double F_inc;
  F_inc = stefan_boltzmann * T_star * T_star * T_star * T_star
    * r_star * r_star / (a_AU*a_AU*AU_in_cm*AU_in_cm);
  return F_inc;
}

double min_incident_angle(double dtot, double R_star)
{
  return (4.0*R_star/(3.0*M_PI*dtot*AU_in_cm));
}

double array_interp(double x, double y, double r0, 
                    double dr, double dphi, double *arr)
{
  // interpolate between points to get value of the array at off-grid spots
   
  int ir, iphi, iphi_p1, ll, ll_p1;
  double r, phi, cr, cphi, val;

  r = cartesian2cylindrical(x,y,&r,&phi);
  
  ir = (int) ((r-r0)/dr);
  iphi = (int) ((phi+M_PI)/dphi);
  iphi_p1 = (iphi+1)%NX;
  

  cr = (r-ir*dr-r0)/dr;
  if(r<r0){
    //in case l.o.s. cuts through the cavity
    cr=0.0;
    ir = 0;
  }
  cphi = (phi-(-M_PI+iphi*dphi))/dphi;

  ll = iphi + ir*Nx;
  ll_p1 = iphi_p1 + ir*Nx;

  // interpolate exponential in r, linear in phi 
  val = (1.0-cphi)*exp((1.0-cr)*log(arr[ll]) + cr*log(arr[ll+Nx]))
       +cphi*exp((1.0-cr)*log(arr[ll_p1]) + cr*log(arr[ll_p1+Nx]));

  return val;
}

double array_interp2(double x, double y, double r0, double *rpos, 
                    double dphi, double *arr, double *arr2, double *val2)
{
  // interpolate between points to get value of the array at off-grid spots
   
  int ir, iphi, iphi_p1, ll, ll_p1, ir_found, dir;
  double r, phi, cr, cphi, val, dr;
  int full_size_y = NY+2*NGHY;

  r = cartesian2cylindrical(x,y,&r,&phi);

  if(r<rpos[0]*5.2) r=rpos[0]*5.2;

  if(((toupper(*SPACING)) == 'L') && ((toupper(*(SPACING+1))) == 'O')){
    dr = (log(rpos[NY+2*NGHY-1]*5.2)-log(rpos[0]*5.2))/(NY+2*NGHY-1);
    ir = (int)((log(r)-log(rpos[0]*5.2))/dr);
    if(ir<0) printf("r=%f YMIN=%f dr=%E\n",r,YMIN*5.2,dr);
    if(ir<full_size_y-1) dr = (rpos[ir+1]-rpos[ir])*5.2;
    else  dr = (rpos[ir]-rpos[ir-1])*5.2;
  }
  else{
    dr = (rpos[1]-rpos[0])*5.2;
    ir = (int) ((r-rpos[0]*5.2)/dr);
  } 
  iphi = (int) ((phi+M_PI)/dphi);
  iphi_p1 = (iphi+1)%NX;
  
  cr = (r-rpos[ir]*5.2)/dr;

  cphi = (phi-(-M_PI+iphi*dphi))/dphi;
  ll = iphi + ir*NX;
  ll_p1 = iphi_p1 + ir*NX;

  // linear interpolation in r and phi
  if(ir<NY+2*NGHY-1){
    val = (1.0-cphi)*((1.0-cr)*(arr[ll]) + cr*(arr[ll+Nx]))
       +cphi*((1.0-cr)*(arr[ll_p1]) + cr*(arr[ll_p1+Nx]));
       
    *val2 = (1.0-cphi)*((1.0-cr)*(arr2[ll]) + cr*(arr2[ll+Nx]))
       +cphi*((1.0-cr)*(arr2[ll_p1]) + cr*(arr2[ll_p1+Nx]));
  }
  else{
    val = (1.0-cphi)*arr[ll]+cphi*arr[ll_p1];
    *val2 = (1.0-cphi)*arr2[ll] +cphi*arr2[ll_p1];
  }

  return val;
}


double integrate_gauss(double x1, double x2, int printint){
  int n;
  int nsteps=10;
  double integ,dtt,tt,a1,a2;
  
  dtt = (x2-x1)/nsteps;
  tt = x1;
  
  a1 = exp(-tt*tt);
  for(n=0;n<nsteps;n++){
    tt+=dtt;
    a2 = exp(-tt*tt);
    integ+=(a1+a2)/2.0;
    a1=a2;
  }
  integ=integ*dtt*1.1283;
  return integ;
}

double cartesian2cylindrical(double rx, double ry, double *r, double *phi)
{
  *r = sqrt(rx*rx + ry*ry);
  if (ry == 0.0){
    if(rx>0) *phi = 0.0;
    else *phi=-M_PI;
  }
  else {
    if (rx == 0.0) {
      if (ry>0) *phi = 0.5*M_PI; 
      else *phi = -0.5*M_PI; 
    }
    else *phi = atan2(ry,rx); //was atan(ry/rx);
  }
  return (*r);
}

double incident_angle(double rmid, double zmid,
		      double dzdr, double dzdphi, double R_star)
{
  /* returns angle of incidence from (0,0) to the point (rmid,zmid)
   * assuming the surface has slopes dzdr and dzdy */
  double mu, dtot;
  double dzdy = dzdphi/rmid;
  dtot = sqrt(rmid*rmid+zmid*zmid);
  mu = (dzdr*rmid-zmid)/(dtot*sqrt(1.0+dzdy*dzdy+dzdr*dzdr));
  if (mu<0.0) mu = 0.0;
  /* add 4*R_star/3*PI*r to account for non-zero flux at normal incidence */
  mu += min_incident_angle(dtot,R_star);
  return mu;

}

double moment_solidangle_patch(point origin, 
			       point p0, point p1, point p2, point p3,
			       vector *unitnormal)
{
  /* calculate the integral of the cosine of the angle with the surface 
   * normal (nu) over the solidangle of the rectangular "patch" 
   * defined by p0, p1, p2, and p3
   *
   * using simplified form for the equation:
   * int nu dOmega = 0.5*[ 
   *   angle0/sin(angle0)*dotproduct(unitnormal,crossproduct(d0,d1))
   * + angle1/sin(angle1)*dotproduct(unitnormal,crossproduct(d1,d2))
   * + angle2/sin(angle2)*dotproduct(unitnormal,crossproduct(d2,d3))
   * + angle3/sin(angle3)*dotproduct(unitnormal,crossproduct(d3,d0))
   * where angle(i) = angle between d(i) and d((i+1)mod 4) 
   * *** THIS ASSUMES THAT ALL 4 POINTS ARE CO-PLANAR ***
   */ 
  vector distvec[SQUARE];
  vector normal;
  point center;
  double angle[SQUARE];
  double cosa[SQUARE];
  double retval;
  int i;

  //center = pointat(0.0,0.0,0.0);

  distvec[0] = unitvector(vectorbtwnpts(origin,p0));
  distvec[1] = unitvector(vectorbtwnpts(origin,p1));
  distvec[2] = unitvector(vectorbtwnpts(origin,p2));
  distvec[3] = unitvector(vectorbtwnpts(origin,p3));

  normal = crossproduct(vectorbtwnpts(p0,p1),vectorbtwnpts(p0,p2));
  *unitnormal = unitvector(normal);
  /* check to see if things oriented properly */
  /** ALLOW BACK SCATTER **
  if (dotproduct(distvec[0],normal)*dotproduct(normal,vectorbtwnpts(p0,center))<=0.0) {
    //    printf("wrong orientation %g %g\n",dotproduct(distvec[0],normal),dotproduct(normal,vectorbtwnpts(p0,center)));
    return 0.0;
  }
  **/

  retval = 0.0;
  for(i=0;i<SQUARE;i++) {
    cosa[i] = dotproduct(distvec[i],distvec[(i+1)%SQUARE]);
    /* check for zeros */
    if (fabs(fabs(cosa[i])-1.0)<1e-15)
      return 0.0;
    angle[i] = acos(cosa[i]);
    retval += 0.5*angle[i]/sin(angle[i])
      *dotproduct(*unitnormal, crossproduct(distvec[i],distvec[(i+1)%SQUARE]) );
  }

  /* BACK SCATTER: return abs val */

  return fabs(retval);
}

double moment_solidangle_strip(double r0, double z0,
			       double rin, double zin,
			       double rout, double zout)
{
  /* returns the integral of nu dOmega over a strip of sky (2D) */
  /* (x0,z0) = point in within the disk
   * (xin,zin) and (xout,zout) define the strip of surface 
   * assume that star is at (0,0,0)
   */
  point origin, p0, p1, p2, p3;
  vector norm;
  double retval;

  /* CYLINDRICAL: approximate strip as a wedge, 
   * rectangular, PI/3 in size */
  origin = pointat(r0,0.0,z0);
  p0 = pointat(rin,rin/sqrt(3.0),zin);
  p1 = pointat(rout,rout/sqrt(3.0),zout);
  p2 = pointat(rout,-rout/sqrt(3.0),zout);
  p3 = pointat(rin,-rin/sqrt(3.0),zin);
  retval = moment_solidangle_patch(origin, p0, p1, p2, p3, &norm);
  return retval;
}

double fit_alpha(int npts, double rvals[], double hvals[])
{
  /* do linear regression to log(rvals) vs log(hvals) to solve for alpha */
  /* see Bevinton Ch 7.4 for calculation of weights */
  double *logr, *logh, *wt;
  int i;
  double norm, alpha;
  double c00, c01, c11, chisq;
  
  logr = (double *) malloc(sizeof(double)*npts);
  logh = (double *) malloc(sizeof(double)*npts);
  wt = (double *) malloc(sizeof(double)*npts);
  for(i=0;i<npts;i++) {
    logr[i] = log(rvals[i]);
    logh[i] = log(hvals[i]);
    wt[i] = hvals[i]*hvals[i];
  }
  gsl_fit_wlinear(logr,1,wt,1,logh,1,npts, 
		 &norm, &alpha, &c00, &c01, &c11, &chisq);

  free(logr);
  free(logh);
  free(wt);

  return alpha;
  
}

int gsl_fit_wlinear (const double *x, const size_t xstride,
                 const double *w, const size_t wstride,
                 const double *y, const size_t ystride,
                 const size_t n,
                 double *c0, double *c1,
                 double *cov_00, double *cov_01, double *cov_11,
                 double *chisq)
{

  /* compute the weighted means and weighted deviations from the means */

  /* wm denotes a "weighted mean", wm(f) = (sum_i w_i f_i) / (sum_i w_i) */

  double W = 0, wm_x = 0, wm_y = 0, wm_dx2 = 0, wm_dxdy = 0;

  size_t i;

  for (i = 0; i < n; i++)
    {
      const double wi = w[i * wstride];

      if (wi > 0)
        {
          W += wi;
          wm_x += (x[i * xstride] - wm_x) * (wi / W);
          wm_y += (y[i * ystride] - wm_y) * (wi / W);
        }
    }

  W = 0;                        /* reset the total weight */

  for (i = 0; i < n; i++)
    {
      const double wi = w[i * wstride];

      if (wi > 0)
        {
          const double dx = x[i * xstride] - wm_x;
          const double dy = y[i * ystride] - wm_y;

          W += wi;
          wm_dx2 += (dx * dx - wm_dx2) * (wi / W);
          wm_dxdy += (dx * dy - wm_dxdy) * (wi / W);
        }
    }

  /* In terms of y = a + b x */

  {
    double d2 = 0;
    double b = wm_dxdy / wm_dx2;
    double a = wm_y - wm_x * b;

    *c0 = a;
    *c1 = b;

    *cov_00 = (1 / W) * (1 + wm_x * wm_x / wm_dx2);
    *cov_11 = 1 / (W * wm_dx2);

    *cov_01 = -wm_x / (W * wm_dx2);

    /* Compute chi^2 = \sum w_i (y_i - (a + b * x_i))^2 */

    for (i = 0; i < n; i++)
      {
        const double wi = w[i * wstride];

        if (wi > 0)
          {
            const double dx = x[i * xstride] - wm_x;
            const double dy = y[i * ystride] - wm_y;
            const double d = dy - b * dx;
            d2 += wi * d * d;
          }
      }

    *chisq = d2;
  }

  return 1;
}
