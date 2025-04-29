#include "fargo3d.h"

const double 
  stefan_boltzmann = 5.6704e-5,
  k_boltzmann = 1.3806e-16,
  meanmolwt = 3.85e-24, /* mu = 2.3 m_H */
  gamma_th = 5.0/3.0,
  AU_in_cm = 1.4960e13,
  year = 3.1557e7,
  grav_const = 6.6743e-8,
  solar_mass = 1.9891e33,
  earth_mass = 5.9736e27;
               

double B_inner_disk_HSL(int j, int j_shift, double *full_surf, double *rpos,
                      double *nu_omega_tot, struct disk_parameters *DP, 
                      struct disk_opacity *opacity, struct var_struct *VS, int cells_in)
{
  //this function calculates the flux from the disk that is interior to
  //the simulation box. Because we don't have detailed information on this
  //region, we assume it to be axisymmetric and dominated by diffuse radiation
  double r0, dr, rlim, redge;
  double B;
  double alpha;
  double *rads,*heights;
  int nfit=NFIT,n,i;
  nfit = 4; 

  redge = rpos[0]*5.2;
  dr =  (rpos[1]-rpos[0])*5.2; 
  rlim = redge/2.0;
  if(rlim<redge-1.5*full_surf[i]) rlim = redge-1.5*full_surf[i];
  if(rlim<redge/10) rlim=redge/10;
  r0 = rpos[j]*5.2;

  //find a power law for the surface of the interior
  heights = (double *)calloc(nfit,sizeof(double));
  rads = (double *)malloc(sizeof(double)*nfit);
  for(n=0;n<nfit;n++) {
    rads[n] = rpos[n]*5.2;
    for(i=0;i<NX;i++){     
      heights[n] += full_surf[i+n*Nx];
    }
    heights[n] *= 1.0/Nx;
  }
  alpha = fit_alpha(nfit,rads,heights);

  double rsurf, zsurf, dzdr;
  double Firr, nu_omega, mu_i;
  double r1, r2, z1, z2;

  //integrate over strips, assuming a power-law for the surface */
  B = 0.0;
  //find surface height following the power law
  r2 = redge;
  z2 = heights[0];
  r1 = redge-dr;
  z1 = z2*pow(r1/r2,alpha);
  
  while(r1>rlim){  
    dzdr = (z2-z1)/(r2-r1);
    rsurf = 0.5*(r1+r2);
    zsurf = 0.5*(z1+z2);      
    mu_i = 4.0e-2; //incident_angle(rsurf,zsurf,dzdr,0.0,DP->R_star);
    nu_omega = moment_solidangle_strip(r0,0.0,r1,z1,r2,z2);
    
    Firr = F_incident(DP->T_star,DP->R_star,sqrt(rsurf*rsurf+zsurf*zsurf));   
    B += nu_omega*B_unpert_diffuse(mu_i,Firr,opacity);
    *nu_omega_tot += nu_omega;
    
    z2 = z1;
    r2 = r1;
    r1 -= dr;
    z1 = z2*pow(r1/r2,alpha);    
  }


  free(heights);
  free(rads);
  return B;
}

double B_outer_disk_HSL(int j, int j_shift, double *full_surf, double *rpos,
                      double *nu_omega_tot, struct disk_parameters *DP, 
                      struct disk_opacity *opacity, struct var_struct *VS)
{
  //this function calculates the flux from the disk that is exterior to
  //the simulation box. Because we don't have detailed information on this
  //region, we assume it to be axisymmetric and dominated by diffuse radiation
  double r0, dr;
  double B;
  double redge, alpha;
  double *rads,*heights;
  int nfit,n,back_j,i;
  nfit = 4; //was 8
  
  back_j = (Ny*CPU_Number+2*NGHY-1);
  r0 = rpos[j]*5.2; 
  //find the power law for the disk exterior
  heights = (double *)calloc(nfit,sizeof(double));
  rads = (double *)malloc(sizeof(double)*nfit);
  for(n=0;n<nfit;n++) {
    for(i=0;i<Nx;i++){
      heights[n] = full_surf[i+(back_j-(nfit-n-1))*Nx];
    }
    heights[n] *= 1.0/Nx;
    rads[n] = rpos[NY+2*NGHY-nfit+n]*5.2;
  }
  alpha = fit_alpha(nfit,rads,heights);
  redge  = rpos[back_j]*5.2;

  double rsurf, zsurf, dzdr;
  double Firr, nu_omega, mu_i, tau;
  double r1, r2, z1, z2;

  dr =  0.03*redge; 

  /*/integrate over strips, assuming a power-law for the surface */
  B = 0.0; 
  // assume a slope consistent with a power-law surface 
  r1 = redge;
  z1 = full_surf[Nx/2+back_j*Nx];//heights[nfit-1];
  r2 = redge+dr;
  z2 = z1*pow(r2/r1,alpha);
  
  while(r2<1.5*redge){  
    dzdr = (z2-z1)/(r2-r1);
    rsurf = 0.5*(r1+r2);
    zsurf = 0.5*(z1+z2);
    
    nu_omega = moment_solidangle_strip(r0,0.0,r1,z1,r2,z2);
    mu_i = incident_angle(rsurf,zsurf,dzdr,0.0,DP->R_star);
    
    Firr = F_incident(DP->T_star,DP->R_star,sqrt(rsurf*rsurf+zsurf*zsurf));   
    B += nu_omega*B_unpert_diffuse(mu_i,Firr,opacity);
    *nu_omega_tot += nu_omega;

    z1 = z2;
    r1 = r2;
    r2 += dr;
    z2 = z1*pow(r2/r1,alpha);    
  }


  free(heights);
  free(rads);
  return B;
}

double B_unpert_diffuse(double mu, double F_irr, 
		    struct disk_opacity *opacity)
{
  //calculate flux from a diffuse source
  double C1, C2;
  double c1p;
  double B;
  double sca_frac, abs_frac, g_param, q;

  sca_frac = opacity->sca_frac;
  abs_frac = opacity->abs_frac;
  g_param = opacity->g_param;
  q = opacity->ratio;

  if (mu<1e-16) {
    return 0.0;
  }

  C1 = -3.0*sca_frac*mu*mu/(1-g_param*g_param*mu*mu);
  C2 = sca_frac*(2.0+3.0*mu)/( g_param*(1.0+2.0*g_param/3.0)
				* (1-g_param*g_param*mu*mu) );
  
  c1p = (1.0+C1)*(2.0+3.0*mu/q) 
    + C2*(2.0+3.0/g_param/q);

  B = abs_frac*F_irr*mu/4.0/M_PI * c1p;
					    
  return B;
}	

void calc_flux_azi5(struct disk_parameters *DP,struct disk_opacity *opacity,struct var_struct *VS,
                   double *full_dens, double *full_H, double *full_surf, float *full_cs,
                   double *flux_azi, double *dens_azi,int opac_i, double *rpos)
{
  //this function calculates the flux at the midplane. A fraction of the
  //midplane points have flux from all of the surface points calculated using
  //the moment method from Jang-Condell 2008. The remaining midplane points
  //have their received flux calculated by interplation from the other points
  
  int i, j, k, ii, jj, ii2, jj2, ll, l2;
  int jmin, jmax, jmin_prev, jmax_prev, jmin2, jmax2;
  int icount, i_shift, imin, imax, i2;
  double r, dr, r0;
  double opac, opacR, nu_omega, flux_sub, nu_omega_sub, flux_sub2, nu_omega_sub2;
  point origin;
  vector vert; 
  double flux_inner, flux_outer, nu_om_out, nu_om_in, tau_min;
  int im2, im1, ip1, ip2;
  double Jmin =  5.67e-5*1e4/4/M_PI;

  
  int size_x = NX/CPU_Number;
  int size_y = NY+2*NGHY;
  int full_size = NX*(NY+2*NGHY);

  i_shift = CPU_Rank*size_x;
  vert.x = 0.0; vert.y=0.0; vert.z=1.0;

  int i_step, j_step, n_i_steps, n_j_steps;
  int *i_steps, *j_steps;
  int no_interp=0;
  i_step = NX/NX_RAD; j_step = (int) NY/NY_RAD;  //64
  if(i_step==0) i_step=1;
  if(j_step==0) j_step=1;
  n_i_steps = size_x/i_step; 
  n_j_steps = size_y/j_step;
  if(n_i_steps<2) n_i_steps=2;
  if(i_step*(n_i_steps-1)<size_x-2) n_i_steps+=1;
  if(j_step*(n_j_steps-1)<size_y-2) n_j_steps+=1;
  i_steps = (int *) calloc(n_i_steps,sizeof(int));
  j_steps = (int *) calloc(n_j_steps,sizeof(int));
  for(i=0;i<n_i_steps-1;i++) i_steps[i] = i*i_step+i_shift;
  i_steps[n_i_steps-1]=size_x-1+i_shift;
  for(i=0;i<n_j_steps-1;i++) j_steps[i] = i*j_step;
  j_steps[n_j_steps-1]=size_y-1;

  //scale free opacity, since the length dl will be scale free
  opac = (opacity->Planck);
  opacR = opacity->Rosseland;

  i=j=k=0;
  origin.z = 0.0;

  double jlim1,jlim2,ilim1,ilim2;

  double tau_return;
  int flux_method;
  double flux_main, lim_min, lim_max, tau_lim, lim_ratio, calc_frac;
  int jp1,jm1,ip3, l3,j3,i3,i4, j3min, j3max, calc_step,l_region,l_region_start;
  int quick1, full1, quick2, full2, skip2, tau2, tot_count, calc_count;
  double dx1,dy1,dz1,dx2,dy2,dz2;
  int side_x1,side_x2,side_y1,side_y2,lenx,leny,region_size;
  double sidex1,sidex2,sidey1,sidey2, midx_ll, midy_ll, dl1, tau_l;
  float Jmin_10K = 1.8e-1; //erg/cm2/s

  quick1=quick2=full1=full2=tau2=skip2=0;
  int point_considered = 0;
  int print_flux_irr=1;

  lim_min = FLXLIMMIN; lim_max = FLXLIMMAX; 
  lim_ratio = log(lim_max/lim_min);

  sidex1 = (NX/NXRADSURF-1)/2.0;
  sidex2 = NX/NXRADSURF/2.0;
  sidey1 = (NY/NYRADSURF-1)/2.0;
  sidey2 = NY/NYRADSURF/2.0;
  side_x1 = (int)(sidex1); side_x2 = (int)(sidex2);
  side_y1 = (int)(sidey1); side_y2 = (int)(sidey2);
  lenx=side_x1+side_x2+1;
  leny=side_y1+side_y2+1;
  region_size = lenx*leny;

 
  point s[4];
  vector v;
  double ph_nx, ph_ny, ph_nz, mag, r_shadow, z_shadow, zph, zph2, max_slope, cos_surf;
  int jmins[n_i_steps][n_j_steps], jmaxs[n_i_steps][n_j_steps];

  int jrange = 2;
  int irange = NX/16;
  int jbuff = 0; //NX/64;

  for(ii2=0; ii2<n_i_steps; ii2++){
    ii = i_steps[ii2];
    jmin_prev = 1;
    jmin2 = 1;
    jmax_prev = size_y-1;
    jmax2 = size_y-1;
    for(jj2=0;jj2<n_j_steps;jj2++){
      jj = j_steps[jj2];
      ll = ii+jj*NX;
      zph = get_photosphere(full_dens[ll], full_H[ll], opacR, 2.0/3, 0.0)*5.2;
      r = rpos[jj]*5.2;

      max_slope = 0.0;
      jmin = 1;
      for(j3=0;j3<jj-jbuff;j3++){
        zph2 = get_photosphere(full_dens[ii+j3*NX], full_H[ii+j3*NX], opacR, 2.0/3, 0.0)*5.2;
        if((zph2-zph)/(jj-j3)>max_slope){
          max_slope = (zph2-zph)/(jj-j3);
        } 
      }
      if(max_slope>0.0){
        jmin = jj- (int) ((full_surf[ll]-zph)/max_slope*cos(atan(full_surf[ll]/r))) ;
        if(jmin<1) jmin=1;
      }
      if(jmin>1 || jmin_prev>1){
        if(jmin>=jmin_prev) jmin2 = (jmin+(jmin_prev+j_step))/2;
        else jmin2 = (3*jmin+jmin_prev+j_step)/4;
      }
      jmin_prev = jmin2;
      jmins[ii2][jj2]=jmin2;  
    }

    max_slope = 0.0;
    jmax = size_y-1;
    for(jj2=n_j_steps-1;jj2>=0;jj2--){
      jj = j_steps[jj2];
      ll = ii+jj*NX;
      zph = get_photosphere(full_dens[ll], full_H[ll], opacR, 2.0/3, 0.0)*5.2;  
      max_slope=0.0;
      jmax=size_y-1;
      for(j3=size_y-1;j3>jj+jbuff;j3--){
        zph2 = get_photosphere(full_dens[ii+j3*NX], full_H[ii+j3*NX], opacR, 2.0/3, 0.0)*5.2;
        if((zph2-zph)/(j3-jj)>max_slope) max_slope = (zph2-zph)/(j3-jj);
      }
      if(max_slope>0.0){
        jmax = jj + (int) ((full_surf[ll]-zph)/max_slope/cos(atan(full_surf[ll]/r))) ;
        if(jmax>size_y-1) jmax=size_y-1;
      }
      if(jmax<size_y-1 || jmax_prev<size_y-1){
        if(jmax<=jmax_prev) jmax2 = (jmax+jmax_prev-j_step)/2;
        else jmax2 = (3*jmax+jmax_prev-j_step)/4;
      }
      jmax_prev = jmax2;
      jmaxs[ii2][jj2]=jmax2;
    }
  }

  double plaw_factor;
  real mu = 2.3;
  real mp = 1.67e-24 ;  //g
  real k_b = 1.381e-16; // erg K-1
  real gamma = 1.667;
  real StfBolt = 5.67e-5;
  real f0,cs0;

  double phi,azi_shadow;
  azi_shadow = 1.0;
#ifdef AZISHADOW
  double shade_angle = M_PI/10;
  double shaded = 0.01;
#endif


  //main loop for calculating the flux
  //loop through midplane elements, but not every one
  for(jj2=n_j_steps-1;jj2>=0;jj2--){
    jj = j_steps[jj2];

    jp1 = jj+1;
    jm1 = jj-1;
    if(jj==0) jm1 = 0;
    if(jj==size_y-1) jp1=jj;
    //calculate the flux from the inner and outer disk beyond the domain.
    //assumed to be axisymmetric, so calculate once at each radius
    nu_om_in=0.0;
    flux_inner=B_inner_disk_HSL(jj,0,full_surf,rpos,&nu_om_in,
                                 DP,opacity,VS,jrange-jj)/M_PI;
      //default size is angular size is Pi/3, but we may not be considering that much
    nu_om_out=0.0;
    flux_outer=B_outer_disk_HSL(jj,0,full_surf,rpos,&nu_om_out,
                                 DP,opacity,VS)/M_PI;                          

    for(ii2=0; ii2<n_i_steps; ii2++){
      ii = i_steps[ii2];
      ll = ii+jj*NX;
      l2 = ii-i_shift + jj*NX/CPU_Number;
      im1 = ii-1;
      ip1 = ii+1;
      if(ii==0) im1 = NX-1;
      if(ii==NX-1) ip1 = 0;
      
      origin.x = VS->midx[ll];
      origin.y = VS->midy[ll]; 
      nu_omega_sub = 0.0;
      flux_sub = 0.0;
      //loop through all surface elements

      i=ii;
      j=jj;
      r = rpos[j]*5.2; 
      dr = (rpos[j+1]-rpos[j])*5.2;
      tau_min = get_tau_min(ii,jj,full_dens,full_H,opac/opacity->ratio,rpos,r0,VS);

#ifdef AZISHADOW
      azi_shadow = 1.0;
      phi = -M_PI+2.0*i/NX*M_PI;
      if(fabs(M_PI/2-fabs(phi))>M_PI/2.0-shade_angle) azi_shadow = shaded;
#endif
      flux_method=0;
      flux_main = B_flux_azi5(VS,full_H,full_dens,i,j,ii,jj,i_shift,ll,
              &nu_omega, NX, origin.x, origin.y, 0.0,
              opacity,r,rpos,r0,&tau_return,&flux_method,tau_min,
              ph_nx, ph_ny, ph_nz, zph)/M_PI*azi_shadow;
        
        //if you want to go full plane parallel
        flux_sub = flux_main;
        nu_omega_sub = nu_omega;

        imin=ii-irange; 
        imax=ii+irange;  
        jmin= jmins[ii2][jj2];
        jmax= jmaxs[ii2][jj2];
        if(jmin<1) jmin=1;
        if(jmax>size_y-1) jmax=size_y-1;
        for(j=jmin;j<jmax;j+=leny){   
          for(i2=imin;i2<imax;i2+=lenx){
            i = (i2+NX)%NX;
            r = rpos[j]*5.2; 
            dr = (rpos[j+1]-rpos[j])*5.2;
#ifdef AZISHADOW
            azi_shadow = 1.0;
            phi = -M_PI+2.0*i/NX*M_PI;
            if(fabs(M_PI/2-fabs(phi))>M_PI/2.0-shade_angle) azi_shadow = shaded;
#endif
            //add the flux from this particular surface cell to the total  
            flux_method=0;       
            flux_sub2 = B_flux_azi5(VS,full_H,full_dens,i,j,ii,jj,i_shift,ll,
                    &nu_omega, NX, origin.x, origin.y, 0.0,
                    opacity,r,rpos,r0,&tau_return,&flux_method,tau_min,
                    ph_nx, ph_ny, ph_nz, zph)/M_PI*azi_shadow;
            nu_omega_sub2 = nu_omega;
            calc_count = 1;

            if(flux_sub2>flux_main){
              flux_main=flux_sub2;
            } 
            if(flux_sub2<lim_min*flux_main){
              flux_sub2*=region_size;
              nu_omega_sub2*=region_size;
            }
            else{
              if(flux_method==1){
                dx1=(VS->midx[ll]-VS->midx[l]);
                dy1=(VS->midy[ll]-VS->midy[l]);
                dz1=VS->midz[l];
                dl1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
                midx_ll = VS->midx[ll];
                midy_ll = VS->midy[ll];      
                tau_l = tau_return;   
              }
              j3min=j-side_y1; if(j3min<1) j3min=1;
              calc_frac = 1.0/(log(flux_sub2/flux_main/lim_min)/lim_ratio);
              calc_step = (int) (1.0/(log(flux_sub2/flux_main/lim_min)/lim_ratio));
              if(calc_step%lenx==0) calc_step++;
              l_region_start = (ii+jj+i+j)%calc_step;
              for(l_region=l_region_start;l_region<region_size;l_region+=calc_step){
                i3 = (l_region%lenx+i-side_x1+NX)%NX;
                j3 = l_region/lenx+j-side_y1;
                if(j3==j && i3==i) i3++;
                if(j3<1) j3=1; if(j3>size_y-1) j3=size_y-1;
                l3 = j3*NX+i3;
                r = rpos[j3]*5.2; 
#ifdef AZISHADOW
                azi_shadow = 1.0;
                phi = -M_PI+2.0*i3/NX*M_PI;
                if(fabs(M_PI/2-fabs(phi))>M_PI/2.0-shade_angle) azi_shadow = shaded;
#endif
                if(flux_method==1){
                  flux_method=5;                  
                  dx2=(midx_ll-VS->midx[l3]);
                  dy2=(midy_ll-VS->midy[l3]);
                  dz2=VS->midz[l3];
                  tau_return = tau_l*sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)/dl1;
                  if(tau_return<0.0) printf("uh oh %d %d,%d->%d,%d tau_l=%.2E dl1=%.2E\n",CPU_Rank,i,j,ii,j,tau_l,dl1);
                }
                flux_sub2 += B_flux_azi5(VS,full_H,full_dens,i3,j3,ii,jj,i_shift,ll,
                  &nu_omega, NX, origin.x, origin.y, 0.0,
                  opacity,r,rpos,r0,&tau_return,&flux_method,tau_min,
                  ph_nx, ph_ny, ph_nz, zph)/M_PI*azi_shadow;
                nu_omega_sub2 += nu_omega;
                calc_count++;
              }// end for loop of additional l points
              flux_sub2 *= 1.0*region_size/calc_count; 
              nu_omega_sub2 *= 1.0*region_size/calc_count;
            }//end else 
            flux_sub += flux_sub2;
            nu_omega_sub += nu_omega_sub2; 
          } //end for i (main loop)
        }  //end for j (main loop)

#ifdef AZISHADOW
        azi_shadow = 1.0;
        phi = -M_PI+2.0*ii/NX*M_PI;
        if(fabs(M_PI/2-fabs(phi))>M_PI/2.0-shade_angle) azi_shadow = shaded;
#endif
        //add the contributions from the inner and outer disk to the total
        if(jmin==1){
          flux_sub += flux_inner*azi_shadow;
          nu_omega_sub += nu_om_in;
        }
        if(jmax==size_y-1){
          flux_sub+= flux_outer*azi_shadow;
          nu_omega_sub+= nu_om_out;
        }

        if(nu_omega_sub>0.0) flux_azi[l2] = flux_sub*M_PI/nu_omega_sub + Jmin_10K; 
        else flux_azi[l2] = Jmin_10K;
        
        /*if(jj2<n_j_steps-1){
          plaw_factor = exp(-pow(rpos[jj2]-YMIN+dr/5.2*NGHY,2)/2/pow((YMAX-YMIN)/40,2));
          l3 = ii-i_shift + j_steps[jj2+1]*NX/CPU_Number;
          cs0  = ASPECTRATIO*pow(rpos[jj]/R0,FLARINGINDEX)*sqrt(G_CGS*MSTAR_CGS/rpos[jj]/R0_CGS);
          f0=StfBolt*pow(cs0*cs0*mu*mp/k_b/GAMMA,4)/M_PI;
          flux_azi[l2] = exp((1-plaw_factor)*log(flux_azi[l2])
                       + plaw_factor*log(f0));
        }*/

    if(isnan(flux_azi[l2]) || isinf(flux_azi[l2]) || flux_azi[l2]<0.0){
      printf("second NaN catch: CPU=%d - %d,%d\n",CPU_Rank,ii,jj);
      printf("           fs=%.2E nuom=%.2E mu=%.2E tau_min=%.2E\n",flux_azi[l2],nu_omega_sub,VS->mu[l],tau_min);
      exit(0);
    }
  
    }
  }

  double w, w1, w2, w3, w4;
  //2 1-D cubic interps of log(flux)
  if(NY_RAD<NY){
    //first interpolate in the radial direction at the fixed aximuths
    jj2=1;
    for(j=1;j<size_y-1;j++){
      if(j==j_steps[jj2]){
        jj2++;
        j++;
      } 
      for(ii2=0;ii2<n_i_steps;ii2++){
        i=i_steps[ii2]-i_shift;
        if(jj2>2 && jj2<n_j_steps-2){
          w1=1.0*(j-j_steps[jj2-1])*(j-j_steps[jj2])*(j-j_steps[jj2+1])
             /(1.0*(j_steps[jj2-2]-j_steps[jj2-1])*(j_steps[jj2-2]-j_steps[jj2])
                *(j_steps[jj2-2]-j_steps[jj2+1]));                   
          w2=1.0*(j-j_steps[jj2-2])*(j-j_steps[jj2])*(j-j_steps[jj2+1])
             /(1.0*(j_steps[jj2-1]-j_steps[jj2-2])*(j_steps[jj2-1]-j_steps[jj2])
                *(j_steps[jj2-1]-j_steps[jj2+1]));                   
          w3=1.0*(j-j_steps[jj2-2])*(j-j_steps[jj2-1])*(j-j_steps[jj2+1])
             /(1.0*(j_steps[jj2]-j_steps[jj2-2])*(j_steps[jj2]-j_steps[jj2-1])
                *(j_steps[jj2]-j_steps[jj2+1]));                   
          w4=1.0*(j-j_steps[jj2-2])*(j-j_steps[jj2-1])*(j-j_steps[jj2])
            /(1.0*(j_steps[jj2+1]-j_steps[jj2-2])*(j_steps[jj2+1]-j_steps[jj2-1])
               *(j_steps[jj2+1]-j_steps[jj2]));
        
          flux_azi[la]=exp(w1*log(flux_azi[j_steps[jj2-2]*size_x+i])+w3*log(flux_azi[j_steps[jj2]*size_x+i])
                     +w2*log(flux_azi[j_steps[jj2-1]*size_x+i])+w4*log(flux_azi[j_steps[jj2+1]*size_x+i]));
        }
        else{
            w3 = 1.0*(j-j_steps[jj2-1]);
            w2 = 1.0*(j_steps[jj2]-j);
            w = w3+w2;
            flux_azi[la]=exp((w2*log(flux_azi[j_steps[jj2-1]*size_x+i])
                               +w3*log(flux_azi[j_steps[jj2]*size_x+i]))/w);
        }
        if(isnan(flux_azi[la]) || isinf(flux_azi[la])){
          printf("error in flux radial interpolation: CPU=%d - %d,%d\n",CPU_Rank,i,j);
          printf("    flux=%.2E f-1=%.2E, w=%.2f    f+1=%.2E, w=%.2f  \n",
             flux_azi[la],flux_azi[j_steps[jj2-1]*size_x+i],w2,flux_azi[j_steps[jj2]*size_x+i],w3);
          if(jj2>2 && jj2<n_j_steps-2)   printf("     and    f-2=%.2E, w=%.2f    f+2=%.2E, w=%.2f  \n",
             flux_azi[j_steps[jj2-2]*size_x+i],w1,flux_azi[j_steps[jj2+1]*size_x+i],w4);
          exit(0);
        }
      }  
    }
  }
  if(NX_RAD<NX){
    //now linearly interpolate in the azimuth at every r
    for(j=0;j<size_y;j++){
      ii2=0;
      for(i=1;i<size_x-1;i++){
        if(i+i_shift>=i_steps[ii2]) ii2++;
        w2 = 1.0*abs(i+i_shift-i_steps[ii2-1]);
        w1 = 1.0*abs(i+i_shift-i_steps[ii2]);
        w = w1+w2;

        flux_azi[la]=exp((w1*log(flux_azi[j*size_x+i_steps[ii2-1]-i_shift])
                       +w2*log(flux_azi[j*size_x+i_steps[ii2]-i_shift]))/w);
      
      } 
    }
  }
  
  free(j_steps);
  free(i_steps);
}


double B_flux_azi5(struct var_struct *VS, double *full_H, double *full_dens,
                     int i, int j, int ii, int jj, int i_shift, int ll,
                     double *nu_omega, int ygrid, double ox, double oy, double oz,
                     struct disk_opacity *opacity, double r, double *rpos, double r0,
                     double *tau_return, int *flux_method, double tau_min,
                     double ph_nx, double ph_ny, double ph_nz, double zph)
{
  //a function to quickly calculate the flux received at one midplane
  //point from one surface point.                
  double nu, nuph, omega, tau, z_exp, H_root2, dens_int, w1, w2, B, B2, B3, tau2;
  int k,m, max_m, m_step=1;
  point midpoint, origin;
  vector d0, normal;
  k=0;
  double d0x,d0y,d0z,d0_mag,d0zph,angle, cutoff, cutoff2,printint=0;
  double normx_mod, normy_mod, normz_mod, norm_mag, nu_mod;
 
  int ir, iphi, iphi_p1, ll2, ll2_p1, im1, ip1, lm1, lp1;
  double r2, phi, cr, cphi, val;
  double erf1,erf0;

  double dtau_prev, dtau, dtau_step, tau1, r1, tau_v, tau_q, tfrac;
  double mu;
  double Tph;

  int full_size; 
  int size_y = (Ny*CPU_Number+2*NGHY);

  double opac = opacity->Planck;
  double opacR = opac/opacity->ratio;
  double g_param = opacity->g_param;
  double I = 0.0;

  d0x = VS->midx[l]-ox;
  d0y = VS->midy[l]-oy;
  d0z = VS->midz[l]-oz;
  d0_mag = sqrt(d0x*d0x+d0y*d0y+d0z*d0z);

  //normx_mod = 0.0;
  //normy_mod = 0.0;
  //normz_mod = 1.0;
  //mod normal vector is vertical, hence just d0z remains in dot product
  nu_mod = d0z/d0_mag;
  nu_mod = nu_mod;
  nu = (VS->normx[l]*d0x+VS->normy[l]*d0y+VS->normz[l]*d0z)/d0_mag; 
  omega = nu*VS->omeg[l]/(d0_mag*d0_mag);//*W;
  *nu_omega = omega*nu_mod;

  d0zph = VS->midz[l]-zph;
  d0_mag = sqrt(d0x*d0x+d0y*d0y+d0zph*d0zph);
  nuph = (ph_nx*d0x+ph_ny*d0y+ph_nz*d0zph)/d0_mag; 
  if(nu<1.0e-2){
    *nu_omega=0.0;
    return 0.0; 
  }

  double c2, c3;
  c2 = c3 = 0.0;
  int run_tau=0;
  double frac = 1000.0;
  double c1p, c2p, c3p, Firr;
  c1p= (VS->c1p[l]); 
  c2p=VS->c2p[l]; c3p=VS->c3p[l];
  mu = (VS->mu[l]); 
  Firr =  (VS->F_irr[l]); 
  if(mu<=0.0){
    return 0.0; 
  }

  
  double tau_v2 = full_dens[l]*opac/2;
  tau = full_dens[l]*opac/2;

  double diminish, zph_l, zph_l_mod, zph_ll, rl, rll, rhol, rholl, z_los;


  *flux_method=4;
  //run_tau=1;
  tau_v = full_dens[ll]*opacR/2;
  if(*flux_method<=1){
    tau_v = (0.99*full_dens[ll]+0.01*full_dens[l])*opacR/2;
    c3 = c3p*exp(-g_param*(tau_v+2/3.0*mu));  //not modified by nu
    c2 = c2p*exp(-tau_v/mu-2/3.0);            //bc tau_v is ~vert
    if((c2+c3)>c1p/frac) run_tau=1;
    tau_v = full_dens[ll]*opacR/2;
  }
  if(run_tau==1){
    //if the optical depth is sufficiently low in the area, then we cannot
    //assume the diffuse radiation dominates and tau should be calculated
    if(i==ii && j==jj){
      tau=tau_v;
      *tau_return = tau;
    }
    else{
      angle = 1.0/(sin(atan(d0z/sqrt(d0x*d0x+d0y*d0y))));  
        tau=dtau_prev=0.0;
        //interpolate rho & H between midplane points along l.o.s.
        max_m = (int)(sqrt((i-ii)*(i-ii)+(j-jj)*(j-jj)));
        if(abs(i-ii)>Nx/2) max_m = (Nx-abs(i-ii))+abs(j-jj);
        max_m = (int)(max_m*nu)+3*sqrt(512.0/NX);
        //calculating los to photosphere
        zph_ll = get_photosphere(full_dens[ll], full_H[ll], opacR, 2.0/3, 0.0);
        z_exp = (VS->midz[l]/5.2-zph_ll);
    
        //calculate the tau cutoff beyond which the non-diffuse flux is negligible
        cutoff = mu*(2.0/3.0-log(c1p/(frac*c2p)))/nu/opac/angle;
        cutoff2=(-mu*2.0/3.0-log(c1p/(frac*c3p))/g_param)/nu/opac/angle;
        if(cutoff2 > cutoff) cutoff=cutoff2;
  
        erf0 = 0.0; 
        for (m=0;m<max_m;m+=m_step){    
          dens_int = array_interp2(ox+d0x*(m+0.5)/max_m,oy+d0y*(m+0.5)/max_m,
                               r0,rpos,2*M_PI/Nx,full_dens,full_H,&H_root2); 
          H_root2 *= 1.4142; 
          
          erf0 = erf((zph_ll+z_exp*m)/H_root2);
          erf1 = erf((zph_ll+z_exp*(m+1))/H_root2);
          dtau = dens_int/2*(erf1-erf0);

          //simply add dtau
          tau+=dtau*opac*angle;

        }
      //}
        if(z_exp<0.0) tau*=-1.0;
    }

    c3 = c3p*exp(-g_param*(nu*tau+2/3.0*mu));
    c2 = c2p*exp(-nu*tau/mu-2/3.0);
    B = VS->F_irr[l]*mu*(c1p+c2+c3);

    *flux_method=1;
    *tau_return=tau;
 
  }
  //this used to be an else if when flux_method=1 was on option I was considering
  if(*flux_method==5){
    tau = *tau_return;
    c3 = c3p*exp(-g_param*(nu*tau+2/3.0*mu));
    c2 = c2p*exp(-nu*tau/mu-2/3.0);
    B = VS->F_irr[l]*mu*(c1p+c2+c3);
  }
  else{
    //in the diffuse radiation case, the flux simplifies nicely
    B = Firr*mu*c1p;
    *flux_method=4;
  }

  return B*nu_mod*omega;
}	  




double get_tau_min(int ii, int jj,double *full_dens,double *full_H,double opac,
                   double *rpos, double r0, struct var_struct *VS){
  int diri, dirj, min_i, min_j, n, l_dir, ncount_min, m, max_m,xtent,ytent,ytent0;
  double tau_min, dens_avg_min, dens_avg, dens0, dens1, tau_avg, dens_min;
  double d0x, d0y, d0z, d0_mag, x0, y0;
  double tau, dtau, z_exp, H_root2, dens_int, erf0, erf1, angle;
  int ncount1 = 10;
  int ncount2 = ncount1;
  int thick_count, count;
  double rad_retained = 0.0;

  xtent = NX/16; 
  ytent = NY/16;

  thick_count = count = 0;

  dens_avg = 0.0;
  dens_min = full_dens[ii+NX*jj];
  ncount1 = 0;
  x0 = VS->midx[ii+NX*jj];
  y0 = VS->midy[ii+NX*jj];
  for(diri=-xtent;diri<xtent+1;diri++){
    ytent0 = (int)(fabs(ytent*sin(acos(1.0*diri/xtent))));
    for(dirj=-ytent0;dirj<ytent0+1;dirj++){
      l_dir = (ii+diri+NX)%NX + (jj+dirj)*NX; 
      if(l_dir>=0 && l_dir<NX*(NY+2*NGHY)){

        d0x = VS->midx[l_dir]-x0;
        d0y = VS->midy[l_dir]-y0;
        d0z = VS->midz[l_dir];
        d0_mag = sqrt(d0x*d0x+d0y*d0y+d0z*d0z);
        angle = 1.0/(sin(atan(d0z/sqrt(d0x*d0x+d0y*d0y)))); 

        dens1 = full_dens[l_dir]*angle;
        dens_avg += dens1;

        if(dens1<dens_min) dens_min=dens1;

        count++;
        if(dens1*opac/2>3.0) thick_count++;
      } 
    }
  }
  dens_avg *= 1.0/count;
  tau_avg = dens_avg*opac/2;
  //dens_avg = (dens_min+dens_avg)/2.0;
  tau_min = dens_min*opac/2;

  if(thick_count<count && thick_count>count/2){
    tau_min = tau_min-log(1.0*(count-thick_count)/thick_count);
  }
  tau_min = full_dens[ii+NX*jj]*opac/2;

  return tau_min;
}


void viscous_accretion_flux(double *flux_sub, struct disk_parameters *DP, 
                            struct disk_opacity *opacity){
  double a_inv_cm,tau_d,T_visc4,F_v, Mdot, Mdot_limit, r,dl,dr;
  double omega, viscJ, opacR;
  int i,j,k,ii;
  int size_x = NX; 
  int size_y = Ny+2*NGHY;
  int block_shift = (Ny*CPU_Rank)*NX;
  double accretion_convert = MSTAR_CGS/R0_CGS*sqrt(G_CGS*MSTAR_CGS/R0_CGS);
  opacR = opacity->Rosseland;

  real *rho[NFLUIDS];
  real *vy[NFLUIDS];
  for (ii=0; ii<NFLUIDS; ii++){
    INPUT(Fluids[ii]->Density);
    rho[ii]  = Fluids[ii]->Density->field_cpu;
    INPUT(Fluids[ii]->Vy);
    vy[ii] = Fluids[ii]->Vy->field_cpu;
  }

  Mdot_limit = DP->M_dot*1.9891e33/year;

  i=j=k=0;
  for(j=0;j<size_y;j++){
    r = ymed(j)*5.2; 
    dl = 30*(r/5.2-YMIN)/(YMAX-YMIN);
    omega = sqrt(G_CGS*MSTAR_CGS/pow(r*1.495e13,3));                       //Keplerian frequency
    a_inv_cm = 1.0/(r*AU_in_cm);
    for(i=0;i<size_x;i++){

      Mdot = fabs(2*M_PI*(r/5.2)*(rho[0][l]*vy[0][l]
             +rho[1][l]*vy[1][l]+rho[2][l]*vy[2][l])*accretion_convert);
      if(Mdot>Mdot_limit) Mdot = Mdot_limit;
      Mdot = Mdot_limit*exp(-dl) + Mdot*(1.0-exp(-dl));
      
      F_v = 0.75*6.67e-8*(DP->M_star*1.9891e33)*Mdot*
        M_1_PI*a_inv_cm*a_inv_cm*a_inv_cm*(1.0 - sqrt(DP->R_star*a_inv_cm));
      T_visc4 = 0.5*F_v/stefan_boltzmann;

      tau_d = opacR*rho[1][l]/2.0;

      viscJ = T_visc4*(0.75*tau_d+0.5)*stefan_boltzmann/M_PI;

      flux_sub[l] =  viscJ;
    }
  }
}
