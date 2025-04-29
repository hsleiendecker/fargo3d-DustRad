//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>
#include <time.h>

void surface_range(struct disk_parameters *DP, 
                 struct disk_opacity *opacity, struct var_struct *VS,
                 double *full_surf, double *full_dens, double *rpos)
{
  int i,j,k;
  double min_ang, angle, angle2, dr, wedges;
  int j_shift,ll, ll2; 
  int size_y = NY+2*NGHY;
  int size_x = NX;

  j_shift = CPU_Rank*(Stride/Nx-2*NGHY);
  i=j=k=0;

  double r0,r,mag,z,omega,dphi,d;
  double mu, sca_frac, abs_frac, g_param, q, p, C1, C2;
  int j_start, j_end, im1, ip1, l_m1_phi, l_p1_phi, l_m1_r, l_p1_r;
  double mu_avg, mu_std, mu_frac;
  double opac, tau_m, opacR;
  double max_mu = 1.0e-2;
  opac = opacity->Planck;
  opacR = opacity->Rosseland;
  point s[4], midpoint;
  vector v1, v2, normal, d0, vert;
  
  vert.x = 0.0; vert.y=0.0; vert.z = 1.0;
  dphi = Xmed(1)-Xmed(0);//LP->phipos[1]-LP->phipos[0];
  
  sca_frac = opacity->sca_frac;
  abs_frac = opacity->abs_frac;
  g_param = opacity->g_param;
  q = opacity->ratio;
  p = opacity->ratio_p;
  
  j_start = 1; 
  j_end   = size_y-1;
  double mu1, mu2,r_perc;
  int nmu1, nmu2;
  mu1 = mu2 = 0.0;
  nmu1=nmu2=0;
  
  for(i=0;i<size_x;i++){
    mu_avg = 0.0;
    for(j=j_start;j<j_end;j++){
      r = rpos[j]*5.2;
      dr = (rpos[j+1]-rpos[j])*5.2;
      r_perc = (r/5.2-YMIN+dr/5.2*NGHY)/(YMAX+2*dr/5.2*NGHY-YMIN);
      z = full_surf[l];
      im1 = (i-1+Nx)%Nx;
      ip1 = (i+1+Nx)%Nx;
      l_m1_phi = im1+j*Nx;
      l_p1_phi = ip1+j*Nx;
      l_m1_r = l-Nx;
      l_p1_r = l+Nx;

      s[0] = cylin2cart(rpos[j+1]*5.2, Xmed(i), full_surf[l_p1_r], 0.0); 
      s[1] = cylin2cart(r, Xmed(im1), full_surf[l_m1_phi], 0.0);
      s[2] = cylin2cart(rpos[j-1]*5.2, Xmed(i), full_surf[l_m1_r], 0.0); 
      s[3] = cylin2cart(r, Xmed(ip1), full_surf[l_p1_phi], 0.0);
      v1 = crossproduct(vectorbtwnpts(s[1],s[3]),vectorbtwnpts(s[0],s[2]));
      mag = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
      VS->normx[l] = v1.x/mag;
      VS->normy[l] = v1.y/mag;
      VS->normz[l] = v1.z/mag;
      normal.x = VS->normx[l];
      normal.y = VS->normy[l];
      normal.z = VS->normz[l];
      
      midpoint = cylin2cart(r, Xmed(i), full_surf[l], 0.0); //LP->phipos[i]
      VS->midx[l] = midpoint.x;
      VS->midy[l] = midpoint.y;
      VS->midz[l] = midpoint.z;
      
      d0 = vectorbtwnpts(midpoint,pointat(0.0,0.0,0.0));
      VS->d0_starx[l] = d0.x;
      VS->d0_stary[l] = d0.y;
      VS->d0_stary[l] = d0.z;

      d = sqrt(r*r+z*z);
      mu = cosangle(normal,d0);// +min_incident_angle(d,DP->R_star);
      if(r_perc<0.02){
        mu1 += mu;
        nmu1 += 1;
      } 
      if(r_perc>0.04 && r_perc<0.06){
        mu2 += mu;
        nmu2 += 1;
      } 
      //damp mu towards the inner edge
      mu = 4.0e-2*exp(-(r_perc+0.01)*20) + cosangle(normal,d0)*(1.0-exp(-(r_perc+0.01)*20));// +min_incident_angle(d,DP->R_star);

      if(mu<0.0) mu=0.0;
      VS->mu[l] = mu;
      mu_avg += mu;
      
      VS->F_irr[l] = F_incident(DP->T_star,DP->R_star,d)
                  *abs_frac/4.0/M_PI; 
                  
      VS->omeg[l] = dphi/2.0*((r+dr/2.0)*(r+dr/2.0)-(r-dr/2.0)*(r-dr/2.0))
                /cosangle(normal,vert); //HSL. this can actually be reduced to a purely radial array
                
      C1 = -3.0*sca_frac*mu*mu
           /(1-g_param*g_param*mu*mu);
      C2 = sca_frac*(2.0+3.0*mu)/(g_param*(1.0+2.0*g_param/3.0)
	      * (1-g_param*g_param*mu*mu) );
  
      //the original semi-infinite slab derivation
      VS->c1p[l] = (1.0+C1)*(2.0+3.0*mu/q)+ C2*(2.0+3.0/g_param/q);
      VS->c2p[l] = (1.0+C1)/mu*(q/p-3.0*mu*mu/q);
      VS->c3p[l] = C2*g_param*(q/p-3.0/(q*g_param*g_param));

      if(VS->c2p[l]<0.0) VS->c2p[l]=0.0;
      if(VS->c3p[l]<0.0) VS->c3p[l]=0.0;
      /*if(VS->c1p[l]<0.0 || VS->c2p[l]<0.0){
        if(C1<-1.0){
          VS->c1p[l] = 0.0;
          VS->c2p[l] = 0.0;
        }
        else if(C2<0.0){
          mu = 0.95/g_param;
          C2 = sca_frac*(2.0+3.0*mu)/(g_param*(1.0+2.0*g_param/3.0)
	             * (1-g_param*g_param*mu*mu) );
          VS->c1p[l] = (1.0+C1)*(2.0+3.0*mu/q)+ C2*(2.0+3.0/g_param/q);
          VS->c3p[l] = C2*g_param*(q/p-3.0/(q*g_param*g_param));
        } 
        else if(q/p<3.0*mu*mu/q){
          VS->c2p[l] = 0.0;
        }
      }
      if(VS->c3p[l]<0.0){
        if(C2<0.0){
          mu = 0.95/g_param;
          C2 = sca_frac*(2.0+3.0*mu)/(g_param*(1.0+2.0*g_param/3.0)
	             * (1-g_param*g_param*mu*mu) );
          VS->c1p[l] = (1.0+C1)*(2.0+3.0*mu/q)+ C2*(2.0+3.0/g_param/q);
          VS->c3p[l] = C2*g_param*(q/p-3.0/(q*g_param*g_param));
        } 
      }*/


      if(VS->c1p[l]<0.0 || VS->c2p[l]<0.0 || VS->c3p[l]<0.0){
        printf("ERROR!!! flux constants are negative at (%d,%d)\n",i,j);
        printf(" Opacities may need to be recalculated or the surface is bad\n\n");
        if(VS->c1p[l]<0.0){
          printf("c1p = %.2E = %.2E+%.2E : C1=%.3E mu=%.3E q(ratio)=%.3E C2=%.3E g_param=%.3E\n",
                VS->c1p[l],(1.0+C1)*(2.0+3.0*mu/q),C2*(2.0+3.0/g_param/q),C1,mu,q,C2,g_param);
          printf("check that C1>-1, Everything else >0\n ");
          printf("(g=sqrt(3*ABSFRACNUM/PLANCK)  q=ratio=PLANCK/ROSSELAND\n");
          if(mu>0.5) printf("mu is pretty big:  surf-1=%f surf+1=%f\n",s[2].z,s[0].z);
        }
        if(VS->c2p[l]<0.0){
          printf("c2p = %.2E = %.2E*%.2E: C1=%E mu=%.3E q(ratio)=%.3E p(ratio_p)=%.3E sca_frac=%.3E g=%.3E\n",
                 VS->c2p[l],(1.0+C1)/mu,(q/p-3.0*mu*mu/q),C1,mu,q,p,sca_frac,g_param);
          printf("check that C1>-1 & q^2/p>3*mu^2\n");
          printf("(q=ratio=PLANCK/ROSSELAND  p=ratio_p=RATIOPNUM/ROSSELAND)\n");
          if(mu>0.5) printf("mu is pretty big:  surf-1=%f surf+1=%f\n",s[2].z,s[0].z);
        }
        if(VS->c3p[l]<0.0){
          printf("c3p=%.2E: C2=%.3E  mu=%.3E g_param=%.3E  q=%.3E  p=%.3E \n",VS->c1p[l],C2,mu,g_param,q,p);
          printf("check that C2>0, g>0, & q^2/p (%.2E) > 3/g^2 (%.2E)\n",q*q/p,3/g_param/g_param);
          printf("(g=sqrt(3*ABSFRACNUM/PLANCK)  q=ratio=PLANCK/ROSSELAND  p=ratiop=RATIOPNUM/ROSSELAND)\n");
          if(mu>0.5) printf("mu is pretty big:  surf-1=%f surf+1=%f\n",s[2].z,s[0].z);
        }
        exit(0);
      }
    }
    //fill in the edges of the radial domain
    VS->mu[i] = VS->mu[i+NX];  VS->mu[i+(size_y-1)*NX] = VS->mu[i+(size_y-2)*NX];
    VS->F_irr[i] = VS->F_irr[i+NX];  VS->F_irr[i+(size_y-1)*NX] = VS->F_irr[i+(size_y-2)*NX];
    VS->c1p[i] = VS->c1p[i+NX];  VS->c1p[i+(size_y-1)*NX] = VS->c1p[i+(size_y-2)*NX];
    VS->c2p[i] = VS->c2p[i+NX];  VS->c2p[i+(size_y-1)*NX] = VS->c2p[i+(size_y-2)*NX];
    VS->c3p[i] = VS->c3p[i+NX];  VS->c3p[i+(size_y-1)*NX] = VS->c3p[i+(size_y-2)*NX];
    VS->omeg[i] = VS->omeg[i+NX];  VS->omeg[i+(size_y-1)*NX] = VS->omeg[i+(size_y-2)*NX];
    VS->normx[i] = VS->normx[i+NX];  VS->normx[i+(size_y-1)*NX] = VS->normx[i+(size_y-2)*NX];
    VS->normy[i] = VS->normy[i+NX];  VS->normy[i+(size_y-1)*NX] = VS->normy[i+(size_y-2)*NX];
    VS->normz[i] = VS->normz[i+NX];  VS->normz[i+(size_y-1)*NX] = VS->normz[i+(size_y-2)*NX];

  }

  //add the midpoints for the far ends which aren't in the above loop
  for(i=0;i<size_x;i++){
    for(j=0;j<j_end+1;j+=j_end){
      r = rpos[j]*5.2;

      midpoint = cylin2cart(r, Xmed(i), full_surf[l], 0.0); //LP->phipos[i]
      VS->midx[l] = midpoint.x;
      VS->midy[l] = midpoint.y;
      VS->midz[l] = midpoint.z;
    }   
  }
  
}

void fix_ghost_cells()
{
  int i,j,k,ii,lpx,lmx,count;
  real *rho[NFLUIDS];
  int size_y = Ny+2*NGHY;
  int jmin,jmax,imin,imax;

  for (ii=0; ii<NFLUIDS; ii++){
    INPUT(Fluids[ii]->Density);
    rho[ii] = Fluids[ii]->Density->field_cpu; 
  }  
  
  i=j=k=0;
  count=0;
  
  jmin=Ny;
  jmax=0;
  imin = Nx;
  imax = 0;
  //if(CPU_Rank==1){
    for(j=0; j<size_y; j++) {
      for(i=0; i<Nx; i++) {
        if(rho[1][l]<0.0){
          //if(j<jmin) jmin=j;
          //if(j>jmax) jmax=j;
          //if(i<imin) imin=i;
          //if(i>imax) imax=i;   
          printf("negative dens %d,%d,%d:%.1E ",CPU_Rank,i,j,rho[1][l]);     
        }
      }
    //} 

  }
  
  for(j=size_y/2; j<size_y; j++) {
    for(i=0; i<Nx; i++) {
      for(ii=0;ii<3;ii++){
        if(rho[ii][l]<0.0){
          count++;
          lpx = j*Nx+(i+1)%Nx;
          lmx = j*Nx+(i-1)%Nx;
          if(rho[ii][lpx]>0.0 && rho[ii][lmx]>0.0){
            rho[ii][l] = sqrt(rho[ii][lpx]*rho[ii][lmx]);
          }
          else if(j!=size_y-1 && rho[ii][l+NX]>0.0){
            rho[ii][l] = sqrt(rho[ii][l+Nx]*rho[ii][l-Nx]);
          }
          else{
            rho[ii][l] = rho[ii][l-Nx]*rho[ii][l-Nx]/rho[ii][l-2*Nx];
          }
        }
      }
    }
  }
  for(j=size_y/2-1; j>-1; j--) {
    for(i=0; i<Nx; i++) {
      for(ii=0;ii<3;ii++){
        if(rho[ii][l]<0.0){
          count++;
          lpx = j*Nx+(i+1+Nx)%Nx;
          lmx = j*Nx+(i-1+Nx)%Nx;
          if(rho[ii][lpx]>0.0 && rho[ii][lmx]>0.0){
            rho[ii][l] = sqrt(rho[ii][lpx]*rho[ii][lmx]);
          }
          else if(j!=0 && rho[ii][l-NX]>0.0){
            rho[ii][l] = sqrt(rho[ii][l+Nx]*rho[ii][l-Nx]);
          }
          else{
            rho[ii][l] = rho[ii][l+Nx]*rho[ii][l+Nx]/rho[ii][l+2*Nx];
          }
        }
      }
    }
  }
  
}

struct disk_opacity new_disk_opacity(int opac_i)
{
  struct disk_opacity opacity;
  //small dust, DSHARP (fig 8 w/ max size 32 micron)
  opacity.Rosseland = ROSSELAND*MSTAR_CGS/MSTAR*R0*R0/R0_CGS/R0_CGS;   //disk temp=100 K    chi_R
  opacity.Planck    = PLANCK*MSTAR_CGS/MSTAR*R0*R0/R0_CGS/R0_CGS;     //star temp=4750 K   chi*_P
  opacity.ratio     = opacity.Planck/opacity.Rosseland;
  opacity.ratio_p   = RATIOPNUM*MSTAR_CGS/MSTAR*R0*R0/R0_CGS/R0_CGS/opacity.Rosseland;     //kappa_P/chi_R
  opacity.abs_frac  = ABSFRACNUM*MSTAR_CGS/MSTAR*R0*R0/R0_CGS/R0_CGS/opacity.Planck;       //kappa_P/
  opacity.sca_frac  = 1-opacity.abs_frac; 
  opacity.g_param   = sqrt(3*opacity.abs_frac);
  
  return opacity;
}

struct disk_parameters new_disk_parameters()
{
  struct disk_parameters dp;
  //int i,j,k;
  dp.M_star  = MSTARRAD;
  dp.R_star  = 6.96e10*RSTARRAD; 
  dp.T_star  = TSTARRAD;
  dp.M_dot   = MDOTRAD;
  dp.a_au=(Ymed(Ny+2*NGHY-1)+Ymed(0))/2.0*5.2; //should be in au
  dp.alpha = ALPHA;
  dp.M_planet = 0.0;
  dp.mu0 = 0.01258; /* diskheight is height of disk surface */
  return dp;
}

struct var_struct new_var_struct(){
  struct var_struct vs;
  int size_x=Nx;
  int size_y=Ny+2*NGHY;
  int full_size_y = Ny*CPU_Number+2*NGHY;
  vs.non_diffuse = (int *)calloc(size_x*size_y, sizeof(int));
  vs.normx = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.normy = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.normz = (double *)calloc(size_x*full_size_y,sizeof(double)); 
  vs.midx = (double *)calloc(size_x*full_size_y,sizeof(double)); 
  vs.midy = (double *)calloc(size_x*full_size_y,sizeof(double)); 
  vs.midz = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.F_irr = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.omeg = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.d0_starx = (double *)calloc(size_x*full_size_y,sizeof(double)); 
  vs.d0_stary = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.d0_starz = (double *)calloc(size_x*full_size_y,sizeof(double)); 
  vs.c1p = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.c2p = (double *)calloc(size_x*full_size_y,sizeof(double)); 
  vs.c3p = (double *)calloc(size_x*full_size_y,sizeof(double));
  vs.mu = (double *)calloc(size_x*full_size_y,sizeof(double));
  return vs;
}

void freeVS(struct var_struct *VS){
  free(VS->non_diffuse);
  free(VS->normx);
  free(VS->normy);
  free(VS->normz);
  free(VS->midx);
  free(VS->midy);
  free(VS->midz);
  free(VS->F_irr);
  free(VS->omeg);
  free(VS->d0_starx);
  free(VS->d0_stary);
  free(VS->d0_starz);
  free(VS->c1p);
  free(VS->c2p);
  free(VS->c3p);
  free(VS->mu);
}

void flip_full_array(double *full_array, int x_azi, int x_full, int y_full)
{
  //function to flip the indices of the arrays used for the full
  //disk by each CPU. Necessary after the redoing the domain decomposition
  int i,j,k,c,shift,l2,i2,j2;
  i=j=k=0;
  int azi_domain = x_azi*y_full;
  double *dummy;
  dummy = (double *)malloc((x_full*y_full)*sizeof(double));

  for(j=0;j<y_full;j++){
    for(i=0;i<x_full;i++){
      shift = i/x_azi;
      i2 = i-shift*x_azi;
      l2 = shift*azi_domain+j*x_azi+i2;
      dummy[l] = full_array[l2];
    }
  }
  for(j=0;j<y_full;j++){
    for(i=0;i<x_full;i++){
      full_array[l] = dummy[l];
    }
  }
  free(dummy);
}

void flip_full_array_float(float *full_array, int x_azi, int x_full, int y_full)
{
  //function to flip the indices of the arrays used for the full
  //disk by each CPU. Necessary after the redoing the domain decomposition
  int i,j,k,c,shift,l2,i2,j2;
  i=j=k=0;
  int azi_domain = x_azi*y_full;
  float *dummy;
  dummy = (float *)malloc((x_full*y_full)*sizeof(float));
  for(j=0;j<y_full;j++){
    for(i=0;i<x_full;i++){
      shift = i/x_azi;
      i2 = i-shift*x_azi;
      l2 = shift*azi_domain+j*x_azi+i2;
      dummy[l] = full_array[l2];
    }
  }
  for(j=0;j<y_full;j++){
    for(i=0;i<x_full;i++){
      full_array[l] = dummy[l];
    }
  }
  free(dummy);
}

void update_temp_height_azi(double *full_flux, 
                           struct disk_parameters *DP, struct disk_opacity *opacity, 
                           real dt, int opac_i)
{
  int i,j,k,ii;
  int size_x = NX; 
  int size_y = Ny+2*NGHY;
  int block_shift = (Ny*CPU_Rank)*NX;

  real mu = 2.3;
  real mp = 1.67e-24 ;  //g
  real k_b = 1.381e-16; // erg K-1
  real gamma = 1.667;
  real StfBolt = 5.67e-5; //erg cm-2 s-2 K-4
  //2.238e26 is (cm/au)^2. now convert from ergs/cm2 to scale free
  real energy_convert = G/G_CGS*pow(MSTAR/MSTAR_CGS,2)*pow(R0_CGS/R0,3);
  real velocity_convert = pow(R0_CGS/R0*G/G_CGS*MSTAR/MSTAR_CGS,0.5);

  double opac;
  opac = opacity->Planck;
  opac*=opacity->ratio;
  double thick;
  
  double avg_de = 0.0;
  int count = 0;
  
  real *rho[NFLUIDS];
  real *H[NFLUIDS];

  for (ii=0; ii<NFLUIDS; ii++){
    rho[ii] = Fluids[ii]->Density->field_cpu;  
    H[ii] = Fluids[ii]->HTherm->field_cpu;
  }

  real cs_cgs, Temp, heating, Temp_min, Temp_limit, cs_limit, cs_10K;
  real d_e, de_min;
  real H_eq, flux_in;
  real omega, r;
#ifdef ISOTHERMAL
  real *cs;         
  cs = Fluids[0]->Energy->field_cpu; 
  real e_i;
#endif
#ifdef ADIABATIC
  real *cs;
  cs = (real *) malloc(sizeof(real *)*NX*(Ny+2*NGHY));
  real *energy_dens[NFLUIDS];
  for (ii=0; ii<NFLUIDS; ii++) energy_dens[ii] = Fluids[ii]->Energy->field_cpu;
  i = j = k = 0;
  for(j=0; j<size_y; j++) {
    for(i=XIM; i<size_x; i++) {
      cs[l] = pow((GAMMA-1.0)*energy_dens[0][l]/rho[0][l],0.5);
    }
  }
#endif

  Temp_min = 10; //K
  cs_10K = sqrt(GAMMA*k_b*Temp_min/mu/mp * G/G_CGS*MSTAR/MSTAR_CGS/R0*R0_CGS);
  real tauv, teff;

  i=j=k=0;
  for(j=0; j<size_y; j++) {
    double avg_de = 0.0;
    int count = 0;
    for(i=XIM; i<size_x; i++) {
      cs_cgs = cs[l]/velocity_convert;
      Temp = cs_cgs*cs_cgs*mu*mp/k_b/GAMMA;
      tauv = opac*rho[opac_i][l]/2;
      teff = 3*tauv/8 + sqrt(3)/4 + 1/4/tauv;
      thick = (1.0-exp(-opac*rho[opac_i][l]));
      flux_in = full_flux[l+block_shift]*M_PI*energy_convert;
      heating= (flux_in - StfBolt*pow(Temp,4))*energy_convert;
      //flux is already heating per area, so multiply by time and get energy/area
      d_e = heating * (dt*SCALE_FREE_T);  

      Temp_limit = pow(flux_in/energy_convert/StfBolt,0.25);
      cs_limit=sqrt(Temp_limit*k_b/mu/mp*GAMMA*G/G_CGS*MSTAR/MSTAR_CGS/R0*R0_CGS);
      

#ifdef ADIABATIC
      energy_dens[0][l] += d_e;
      if(energy_dens[0][l]<0.0) energy_dens[0][l] = 0.0;
      cs[l] = pow((GAMMA-1.0)*energy_dens[0][l]/rho[0][l],0.5);
      if((d_e<0.0 && cs[l]<cs_limit) || (d_e>0.0 && cs[l]>cs_limit)){
        cs[l] = cs_limit;
        energy_dens[0][l] = cs[l]*cs[l]*rho[0][l]/(GAMMA-1.0);
      }
      if(cs[l]<cs_10K){
        cs[l]=cs_10K;
        energy_dens[0][l] = cs[l]*cs[l]*rho[0][l]/(GAMMA-1.0);
      } 
#endif
#ifdef ISOTHERMAL
      e_i = cs[l]*cs[l]*rho[0][l]/(GAMMA-1.0);
      e_i += d_e;
      if(e_i<0.0) e_i=0.0;
      cs[l] = pow((GAMMA-1.0)*(e_i)/rho[0][l],0.5);
      //} 
      if(d_e<0.0 && cs[l]<cs_limit){
        //printf("(%d,%d) cs=%.3E -> %.3E (%.3E)\n",i,j,cs[l],cs_min,cs_min/sqrt(G/G_CGS*MSTAR/MSTAR_CGS/R0*R0_CGS));
        cs[l] = cs_limit;
      } 
      else if(d_e>0.0 && cs[l]>cs_limit){
        //printf("(%d,%d) cs=%.3E -> %.3E (%.3E)\n",i,j,cs[l],cs_min,cs_min/sqrt(G/G_CGS*MSTAR/MSTAR_CGS/R0*R0_CGS));
        cs[l] = cs_limit;
      } 
      if(cs[l]<cs_10K) cs[l]=cs_10K;
#endif
      //if(CPU_Rank==0 && i==0 && j==0) printf(" *in rad equilib* ");
      //if(PhysicalTime<2.0) cs[l] = cs_limit;

      r = Ymed(j);
      omega = sqrt(G*MSTAR*MSTARRAD/(r*r*r));
      H_eq = cs[l]/omega; 
      //depends on dt<sound-crossing time
      //using sound-crossing time (1/omega) instead of dynamic (2pi/omega)
      H[0][l] += (H_eq-H[0][l]) * dt * omega;  

#ifdef ADIABATIC
      energy_dens[1][l] = 0.0; //cs[l]*cs[l]*rho[1][l]/(GAMMA-1.0);
      energy_dens[2][l] = 0.0; //cs[l]*cs[l]*rho[2][l]/(GAMMA-1.0);

      if(isnan(energy_dens[0][l]) || isinf(energy_dens[0][l]) || energy_dens[0][l]<=0.0){
        printf("(%d,%d,%d): cs=%.3E e=%.3E f=%.3E rho=%.3E thick=%.2E heat=%.3E  \n  d_e=%.2E  cs_lim=%.3E T_lim=%.3E Temp=%.2E\n",
        CPU_Rank,i,j, cs[l],energy_dens[0][l],flux_in/energy_convert,rho[0][l],thick,heating,
        d_e,cs_limit/sqrt(G/G_CGS*MSTAR/MSTAR_CGS/R0*R0_CGS),Temp_limit,Temp);
        exit(0);
      } 
#endif
    }
  }
#ifdef ADIABATIC
  free(cs);
#endif

}

void thermal_relaxation(real dt)
{
  int i,j,k,ii;
  int size_x = NX; 
  int size_y = Ny+2*NGHY;
  i=j=k=0;

  real mu = 2.3;
  real mp = 1.67e-24;  //g
  real k_b = 1.381e-16; // erg K-1
  real gamma = 1.667;
  real StfBolt = 5.67e-5; //erg cm-2 s-2 K-4
  real energy_convert = G/G_CGS*pow(MSTAR/MSTAR_CGS,2)*pow(R0_CGS/R0,3);
  real velocity_convert = pow(R0_CGS/R0*G/G_CGS*MSTAR/MSTAR_CGS,0.5);

  real e0, cs0, rho0, beta, r, omega, cs, d_e, heating, Temp, Temp0;
  beta = BETACOOLING;
  real *energy_dens[NFLUIDS];
  real *rho[NFLUIDS];
  for (ii=0; ii<NFLUIDS; ii++){ 
    INPUT(Fluids[ii]->Energy);
    OUTPUT(Fluids[ii]->Energy);
    energy_dens[ii] = Fluids[ii]->Energy->field_cpu;
    INPUT(Fluids[ii]->Density);
    rho[ii] = Fluids[ii]->Density->field_cpu;
  }

  double azi_shadow = 1.0;
#ifdef AZISHADOW
  real phi;
  double shade_angle = M_PI/15;
  double shaded = 0.3;
#endif

  //simple thermal relaxation
  for(j=0; j<size_y; j++) {
    r     = Ymed(j);
    omega = sqrt(G*MSTAR/r/r/r);
    cs0 = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;
    //rho0  = SIGMA0*pow(r/R0,-SIGMASLOPE);
    for(i=XIM; i<size_x; i++) {
#ifdef AZISHADOW
      azi_shadow = 1.0;
      phi = Xmed(i);
      if(fabs(M_PI/2-fabs(phi))>M_PI/2.0-shade_angle) azi_shadow = shaded;
#endif
      e0 = cs0*cs0*rho[0][l]/(GAMMA-1.0)*azi_shadow;
      d_e = (energy_dens[0][l]-e0)*omega/(beta*sqrt(YMAX/r)) * dt;
      if(fabs(d_e)>fabs((energy_dens[0][l]-e0))) d_e = (energy_dens[0][l]-e0);
      energy_dens[0][l] -= d_e;
      //energy_dens[0][l] = e0;
    }
  }

}


void _radtransfer_azi_cpu(real dt){
  //main function for the radiative tranfer module
  //overall, the module recalculates the disk's surface, then uses that
  //to calculate the flux from the surface elements to find the heating
  //from stellar radiation. It then updates the energy and scale height

  //////////////////////////////////////////
  // 1. SETUP & REDO DOMAIN DECOMPOSITION //
  //////////////////////////////////////////
  struct disk_opacity opacity;
  struct disk_parameters DP;
  struct var_struct VS;

  int i,j,k,ii,it;
  int size_y = Ny+2*NGHY; 
  int size_x = Nx;
  int full_size_y = Ny*CPU_Number+2*NGHY;
  real r, omega;
  
  //which dust species to use for optical depth
  int opac_i = 1;
  int print_extra = 1;
  
  real *rho[NFLUIDS];
  real *H[NFLUIDS];
  real *e[NFLUIDS];
  for (ii=0; ii<NFLUIDS; ii++){
    INPUT(Fluids[ii]->Density);
    rho[ii]  = Fluids[ii]->Density->field_cpu;
    INPUT(Fluids[ii]->HTherm);
    H[ii]  = Fluids[ii]->HTherm->field_cpu;
    INPUT(Fluids[ii]->Energy);
    e[ii]  = Fluids[ii]->Energy->field_cpu;
  }

  for(j=0;j<size_y;j++){
    for(i=0;i<size_x;i++){
      if(H[0][l]>1e8 || H[0][l]<1.0e-13){
        printf("Correcting H: CPU=%d - H[%d]=%E ->",CPU_Rank,l,H[0][l]);
        if(l<size_x*size_y-1){
          printf("     set H[%d]=H[%d]=%E\n",l,l+1,H[0][l+1]);
          H[0][l]=H[0][l+1];  
        }     
        else{
          printf("     set H_azi[%d]=H[%d]=%E\n",l,l-1,H[0][l-1]);
          H[0][l]=H[0][l-1];  
        } 
      }
    }
  }

  fix_ghost_cells();

  opacity = new_disk_opacity(opac_i);
  DP = new_disk_parameters();

  //Determine the average density at inner edge before redoing the domain
  //decomposition. Needed for the surface calculation.
  double avg_rho0, avg_alpha_rho;
  int di;
  if(CPU_Rank==0){
    di = 10;
    if(NY/CPU_Number<10) di=NY/CPU_Number;
    avg_rho0=0.0;
    avg_alpha_rho=0.0;
    for(i=0;i<NX;i++){
      avg_rho0+=rho[opac_i][i+NGHY*NX];
      avg_alpha_rho = (log(rho[opac_i][i+(NGHY+di)*NX])-log(rho[opac_i][i+(NGHY)*NX]))
                      /(log(Ymed(NGHY+di))-log(Ymed(NGHY)));
    } 
    avg_rho0 *= 1.0/NX;
    avg_alpha_rho*= 1.0/NX;
    //avg with the original slope, but SIGMASLOPE has an inverse sign
    avg_alpha_rho = 0.5*avg_alpha_rho - 0.5*SIGMASLOPE;
  }
  MPI_Bcast(&avg_rho0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&avg_alpha_rho,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  //setup for communicating the parameters from each cpu
  //into new radially sliced domains
  int *recvcounts2;
  int Nx2 = NX/CPU_Number;
  int block_size = NY*Nx2/CPU_Number;
  int ghost_block = Nx2*NGHY;
  int domain_size = full_size_y*NX/CPU_Number;
  recvcounts2 = (int *)malloc(CPU_Number*sizeof(int));
  recvcounts2[0] = block_size+ghost_block;
  for(i=1;i<CPU_Number;i++) recvcounts2[i] = block_size;
  recvcounts2[i-1] += ghost_block;

  int *displs2;
  displs2 = (int *)malloc(CPU_Number*sizeof(int));
  displs2[0] = 0;
  for(i=1;i<CPU_Number;i++) displs2[i] = displs2[i-1]+recvcounts2[i-1];

  if(CPU_Rank==0) block_size+=ghost_block;
  if(CPU_Rank==CPU_Number-1) block_size+=ghost_block;
  double *sub_dens2;
  sub_dens2 = (double *)malloc((block_size)*sizeof(double));
  double *dens_azi;
  dens_azi = (double *)malloc((domain_size)*sizeof(double));
  double *sub_H2;
  sub_H2 = (double *)malloc((block_size)*sizeof(double));
  double *H_azi;
  H_azi = (double *)malloc((domain_size)*sizeof(double));
  double *surf_azi;
  surf_azi = (double *)malloc((domain_size)*sizeof(double));
  float *sub_cs;
  sub_cs = (float *)malloc((block_size)*sizeof(float));
  float *cs_azi;
  cs_azi = (float *)malloc((domain_size)*sizeof(float));

  double *flux_radial;
  flux_radial =  (double *)calloc(Nx*(Ny+2*NGHY), sizeof(double));

  double *flux_azi;
  flux_azi = (double *)malloc((domain_size)*sizeof(double));
  double *sub_flux;
  sub_flux =  (double *)calloc(block_size, sizeof(double));
  

  int nc,nx,ny,ng;
  if(CPU_Rank==0) ng=0;
  else ng=NGHY;
  for(nc=0;nc<CPU_Number;nc++){
    for(nx=0;nx<Nx2;nx++){
      for(ny=0;ny<block_size/Nx2; ny++){
        sub_dens2[nx+ny*Nx2] = rho[opac_i][nx+nc*Nx2+(ny+ng)*Nx];
        sub_H2[nx+ny*Nx2] = H[0][nx+nc*Nx2+(ny+ng)*Nx];
        sub_flux[nx+ny*Nx2] = flux_radial[nx+nc*Nx2+(ny+ng)*Nx];
#ifdef ISOTHERMAL
        sub_cs[nx+ny*Nx2] = e[0][nx+nc*Nx2+(ny+ng)*Nx];
#endif
#ifdef ADIABATIC
        sub_cs[nx+ny*Nx2] = sqrt((GAMMA-1.0)*e[0][nx+nc*Nx2+(ny+ng)*Nx]/rho[0][nx+nc*Nx2+(ny+ng)*Nx]);
#endif
        if(sub_dens2[nx+ny*Nx2]<0.0){
          printf("-");
          sub_dens2[nx+ny*Nx2]=sub_dens2[nx+ny*Nx2-1];
        }
      }
    }
    MPI_Gatherv(sub_dens2, block_size, MPI_DOUBLE, dens_azi, recvcounts2,
            displs2, MPI_DOUBLE, nc, MPI_COMM_WORLD); 
    MPI_Gatherv(sub_H2, block_size, MPI_DOUBLE, H_azi, recvcounts2,
            displs2, MPI_DOUBLE, nc, MPI_COMM_WORLD);
    MPI_Gatherv(sub_cs, block_size, MPI_FLOAT, cs_azi, recvcounts2,
            displs2, MPI_FLOAT, nc, MPI_COMM_WORLD);
    MPI_Gatherv(sub_flux, block_size, MPI_DOUBLE, flux_azi, recvcounts2,
            displs2, MPI_DOUBLE, nc, MPI_COMM_WORLD);
  } 

  free(displs2);
  free(sub_H2);
  free(recvcounts2);
  free(sub_dens2);
  free(sub_cs);
  free(sub_flux);
  free(flux_radial);

  double *rpos;
  rpos = (double *)malloc(full_size_y*sizeof(double));
  double dr;
  if (((toupper(*SPACING)) == 'L') && ((toupper(*(SPACING+1))) == 'O')){
    dr = (log(YMAX)-log(YMIN))/NY;
    for(j=0;j<NY+2*NGHY;j++){
      rpos[j] = exp(log(YMIN) + dr*(j-NGHY));
    }
  } 
  else{
    dr = (YMAX-YMIN)/NY;
    for(j=0;j<NY+2*NGHY;j++){
      rpos[j] = YMIN + dr*(j-NGHY+0.5);
    }
  }

  ////////////////////////////
  // 2. SURFACE CALCULATION //
  ////////////////////////////
  calc_surf_twotemp(&DP,&opacity, dens_azi,H_azi,surf_azi, opac_i,Nx2,full_size_y,avg_rho0,avg_alpha_rho,rpos); //

  /////////////////////////////////////////////////////////////////
  // 3. COMMUNICATE SURFACE, DENSITY, & SCALE HEIGHT TO EACH CPU //
  /////////////////////////////////////////////////////////////////

  int *displs3;
  displs3 = (int *)malloc(CPU_Number*sizeof(int));
  for(i=0;i<CPU_Number;i++) displs3[i] = i*domain_size;
  int *recvcounts3;
  recvcounts3 = (int *)malloc(CPU_Number*sizeof(int));
  for(i=0;i<CPU_Number;i++) recvcounts3[i] = domain_size;
  double *full_surf2;
  full_surf2 = (double *)malloc((size_x*full_size_y)*sizeof(double));
  double *full_dens2;
  full_dens2 = (double *)malloc((size_x*full_size_y)*sizeof(double));
  double *full_H2;
  full_H2 = (double *)malloc((size_x*full_size_y)*sizeof(double));
  float *full_cs;
  full_cs = (float *)malloc((size_x*full_size_y)*sizeof(float));
  

    MPI_Gatherv(surf_azi, domain_size, MPI_DOUBLE, full_surf2, recvcounts3,
                displs3,MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    MPI_Gatherv(dens_azi, domain_size, MPI_DOUBLE, full_dens2, recvcounts3,
                displs3,MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    MPI_Gatherv(H_azi, domain_size, MPI_DOUBLE, full_H2, recvcounts3,
                displs3,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(cs_azi, domain_size, MPI_FLOAT, full_cs, recvcounts3,
                displs3,MPI_FLOAT, 0, MPI_COMM_WORLD);
    

  MPI_Bcast(full_surf2,size_x*full_size_y,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(full_dens2,size_x*full_size_y,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(full_H2,size_x*full_size_y,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(full_cs,size_x*full_size_y,MPI_FLOAT,0,MPI_COMM_WORLD);

  flip_full_array(full_surf2,Nx2,size_x,full_size_y);
  flip_full_array(full_H2,Nx2,size_x,full_size_y);
  flip_full_array(full_dens2,Nx2,size_x,full_size_y);
  flip_full_array_float(full_cs,Nx2,size_x,full_size_y);

  int out_num;
  if(CPU_Rank==0 && print_extra==1){ 
    out_num = (int)((PhysicalTime)/DT)/NINTERM;
    if(PhysicalTime-out_num*NINTERM*DT<2*dt){
      printf(" print surf(%d) ",out_num);
      FILE *fp2;
      char str[150];
      sprintf(str,"%sgathered_surf_full%d.txt",OUTPUTDIR,out_num);
      fp2 = fopen(str, "w+");
      i=j=k=0;
      for(i=0;i<NX;i++){
        for(j=NGHY;j<NY+NGHY;j++){
          fprintf(fp2, "%.5f ",full_surf2[l]);
        }
        fprintf(fp2,"\n");
      }
      fclose(fp2);
    }
  }
  
  VS = new_var_struct();

  ////////////////////////////////////
  // 4. CALCULATE FLUX FROM SURFACE //
  ////////////////////////////////////

  //calculate in advance the values that are used repeatedly by surface cells
  surface_range(&DP,&opacity,&VS,full_surf2,full_dens2,rpos);    

  //calculate the flux from stellar irradation and accretion heating
  calc_flux_azi5(&DP,&opacity,&VS,full_dens2,full_H2,full_surf2,full_cs,
                 flux_azi,dens_azi,opac_i,rpos);

  //////////////////////////////////
  // 5. APPLY THE HEATING/COOLING //
  //////////////////////////////////

  //send the info for the flux of the full disk out to each CPU
  double *full_flux;
  full_flux = (double *)malloc((size_x*full_size_y)*sizeof(double));
  MPI_Gatherv(flux_azi, domain_size, MPI_DOUBLE, full_flux, recvcounts3,
                displs3,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(full_flux,size_x*full_size_y,MPI_DOUBLE,0,MPI_COMM_WORLD);
  flip_full_array(full_flux,Nx2,size_x,full_size_y); 

  if(CPU_Rank==0 && print_extra==1 && PhysicalTime-out_num*NINTERM*DT<2*dt){ 
    printf(" print flux (%d) ",out_num);
    FILE *fp3;
    char str3[100];
    sprintf(str3,"%sgathered_flux%d.txt",OUTPUTDIR,out_num);
    fp3 = fopen(str3, "w+");
    i=j=k=0;
    for(i=0;i<NX;i++){
      for(j=NGHY;j<NY+NGHY;j++){
        fprintf(fp3, "%.3E ",full_flux[l]);
      }
      fprintf(fp3,"\n");
    }
    fclose(fp3);
  }

  //update the energy and scale height (in the original domain decomposition)
  update_temp_height_azi(full_flux,&DP,&opacity,dt,opac_i);

  free(full_flux);
  free(full_surf2);
  free(full_H2);
  free(full_dens2);
  free(full_cs);
  free(cs_azi);
  free(rpos);

  free(flux_azi);
  free(recvcounts3);
  free(displs3);
  free(dens_azi);
  free(H_azi);
  free(surf_azi);
  freeVS(&VS);

}

void RadTransfer(real dt){
  FARGO_SAFE(_radtransfer_azi(dt));
}


void DustEnergyCorrection(){
  int i, j, k, ii;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;

  real *energy_dens[NFLUIDS];
  for (ii=0; ii<NFLUIDS; ii++) energy_dens[ii] = Fluids[ii]->Energy->field_cpu;

#ifdef Z
  for(k=0; k<size_z; k++) {
#endif
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=XIM; i<size_x; i++){
#endif
        energy_dens[1][l]=0.0;
        energy_dens[2][l]=0.0;
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif

}
