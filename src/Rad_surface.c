#include "fargo3d.h"

void surface_smoothing(double *surf_azi){
  if (((toupper(*SPACING)) == 'L') && ((toupper(*(SPACING+1))) == 'O')){
    printf("surface smoothing not yet tested for logarithmic scalings!\n");
  } 
  int m,i,j,l2;
  int size_x = NX/CPU_Number;
  int size_y = NY+2*NGHY;
  int i_shift = CPU_Rank*size_x;
  int box_max = 5;
  int box_len;

  float surf_avg;
  float *smoothed;
  smoothed = (float *)malloc((size_y)*sizeof(float));

  for(i=0;i<size_x;i++){
    for(j=0;j<size_y;j++){
      surf_avg = 0.0;
      if(j>=box_max && j<size_y-box_max) box_len = box_max;
      else if(j<box_max) box_len = j;
      else box_len = size_y-1-j;
      for(m=-box_len;m<box_len+1;m++){
        surf_avg += surf_azi[i+(j+m)*size_x];
      }
      smoothed[j] = surf_avg/(2*box_len+1.0);
    }
    for(j=0;j<size_y;j++){
      surf_azi[i+j*size_x] = smoothed[j];
    }
  }

  free(smoothed);
}



void calc_surf_twotemp(struct disk_parameters *DP, struct disk_opacity *opacity,
                   double *dens_azi, double *H_azi, double *surf_azi,int opac_i,
                   int size_x,int size_y,double avg_rho_in, double avg_alpha_rho, double *rpos)
{
  //A function to calculate the surface of the disk. The optical depth is 
  //calculated along fixed rays and the exact height of the tau=2/3 surface 
  //for the starlight is found by interpolation at the radial grid points
  double rho0, rho1, rho2, tau, dtheta;
  double theta, theta_prev, tan_theta; 
  double opac, dl, Si_subl_r;
  int ii,i,i2,j,k,j2;
  int j_shift=(Stride*CPU_Rank/Nx-2*NGHY*CPU_Rank);
  int surf_found, trapped_count,j_start,j_end,n;
  double new_max_theta=0.0;
  double A,S,F, c1, c2, c3, c4, g1, g2, t_prev, s_prev;
  double dr;  //scale-free
  double r0;  //scale-free (+.5 bc Ymed was used)
  double zph1, zph2, csph, rhoph1, rhoph2, rho2o, Hph1, Hph2, Tph, z, theta_ph;
  double sigma_alt, reduction;
  double energy_convert = G/G_CGS*pow(MSTAR/MSTAR_CGS,2)*pow(R0_CGS/R0,3);
  double velocity_convert = pow(R0_CGS/R0*G/G_CGS*MSTAR/MSTAR_CGS,0.5);
  real mu = 2.3;
  real mp = 1.67e-24 ;  //g
  real k_b = 1.381e-16; // erg K-1
  real gamma = 1.667;
  real StfBolt = 5.67e-5; //erg cm-2 s-2 K-4
  double abs_frac = opacity->abs_frac;
  double zph_depth = 0.1; 
  double surf_depth = 2.0/3;
  double interior_factor=2.5;
  
  //parameters used to extrapolate the density inside the grid inner edge
  S = SIGMASLOPE;
  F = FLARINGINDEX;
  //rough estimate of Silica sublimation radius
  Si_subl_r = DP->T_star*DP->T_star/1250.0/1250.0*DP->R_star*R0/R0_CGS;
  if(YMIN<Si_subl_r){
    printf("Error: inner edge (%f) < estimate for Si sublimation (%f)\n",
           YMIN,Si_subl_r);
    prs_exit(EXIT_FAILURE);
  }
  double dr0;
  double sigma1, sigma2, H1, H2;
  int just_integrate=1;
  int use_gamma=0;
  int jj;
  int nstep0=50;
  dr0 = (rpos[NGHY]-Si_subl_r)/nstep0;


  real omega, r;

  int nray= 60;//y_tot/12;
  double theta_factor;
  double *theta_arr;
  double *tau_arr;
  theta_arr = (real *) malloc(sizeof(double *)*nray);
  tau_arr = (real *) malloc(sizeof(double *)*nray);
  double *theta_arr0;
  double *tau_arr0;
  theta_arr0 = (real *) malloc(sizeof(double *)*nray);
  tau_arr0 = (real *) malloc(sizeof(double *)*nray);

  j_start = NGHY;
  j_end = size_y-NGHY;
  i=k=0;

  //opacity already scale free
  opac = (opacity->Planck);

  i=0;
    ////////////////////////////////////////////////
    // step 1: find the inner edge of the surface //
    ////////////////////////////////////////////////

    //Theta values are selecting using a bisection method to eventually
    //find the theta value for which tau=2/3.
    j=NGHY;
    theta = 0.2;
    theta_prev = 0.0;
    surf_found = 0;
    trapped_count = 0;
    r = rpos[j]; 
    while(surf_found < 1){     
      if(use_gamma==1){
        //the analytic integration of optical depth in the disk interior
        //as a function of r and z results in Gamma functions.
        //These gamma functions are calculated numerically here.
        A = H_azi[la]/r; //ASPECTRATIO;
        c1 = opac*avg_rho_in*pow(r,S-F)/(2.5066*A);
        c2 = S+F+1;
        c3 = 2*F;
        c4 = tan(theta)*tan(theta)*pow(r,c3)/(2*A*A);        
        g1 = calc_upper_incomp_gamma((c2-1)/c3,c4*pow(r,-c3));
        g2 = calc_upper_incomp_gamma((c2-1)/c3,c4*pow(Si_subl_r,-c3));        
        tau = c1*pow(c4,(1-c2)/c3)/c3*(g1-g2);
      }
      else if(just_integrate==1){
        r0 = Si_subl_r;
        tau=0;
        //sigma1 = avg_rho_in*pow(r0/R0,avg_alpha_rho)*(1-r0/rpos[NGHY]) +
          //        dens_azi[la]*pow(r0/rpos[NGHY],avg_alpha_rho)*(r0/rpos[NGHY]);
        sigma1 = avg_rho_in*pow(r0/R0,avg_alpha_rho);
        H1 = ASPECTRATIO*pow(r0/R0,FLARINGINDEX)*r0*(1-r0/rpos[NGHY]) +
                 H_azi[la]*pow(r0/rpos[NGHY],FLARINGINDEX+1)*(r0/rpos[NGHY]);
        zph1 = get_photosphere(sigma1,H1,opacity->Rosseland,zph_depth,0.0) * sqrt(rpos[NGHY]/r0);
        theta_ph=atan(zph1/r0);
        Tph = TSTARRAD*sqrt(RSTARRAD*4.65e-3/r0/cos(theta_ph)/5.2/2)
             /pow(RATIOPNUM/ABSFRACNUM,0.25);      // K
        csph = sqrt(Tph*GAMMA*k_b/mu/mp); // cm/s
        csph *= velocity_convert; // scale free
        Hph1 = csph/sqrt(G*MSTAR*MSTARRAD/pow(r0,3));  //scale free
        rhoph1 = sigma1/(2.5066*H1)*
                exp(-zph1*zph1/2*(1/H1/H1-1/Hph1/Hph1)); 
        sigma_alt = sigma1/2*erf(zph1/1.414/H1) + rhoph1*sqrt(3.1415/2)*Hph1*(1-erf(zph1/1.414/Hph1));
        reduction = sigma1/2/(sigma_alt);
        if(r0*theta>zph1){
          rhoph1*=reduction;
          rho2 = rhoph1*exp(-pow((r0*tan(theta))/(1.414*Hph1),2.0));
        } 
        else{
          sigma1*=reduction;
          rho2 = sigma1/(2.507*H1)*exp(-pow((r0*tan(theta))/(1.414*H1),2.0));
        }
        for(jj=0;jj<nstep0;jj++){
          rho1 = rho2;
          //sigma2 = avg_rho_in*pow((r0+dr0)/R0,avg_alpha_rho)*(1-(r0+dr0)/rpos[NGHY]) +
            //       dens_azi[la]*pow((r0+dr0)/rpos[NGHY],avg_alpha_rho)*((r0+dr0)/rpos[NGHY]);
          sigma2 = avg_rho_in*pow((r0+dr0)/R0,avg_alpha_rho);
          H2 = ASPECTRATIO*pow((r0+dr0)/R0,FLARINGINDEX)*(r0+dr0)*(1-(r0+dr0)/rpos[NGHY]) +
                   H_azi[la]*pow((r0+dr0)/rpos[NGHY],FLARINGINDEX+1)*((r0+dr0)/rpos[NGHY]);
          zph2 = get_photosphere(sigma2,H2,opacity->Rosseland,zph_depth,0.0)* sqrt(rpos[NGHY]/(r0+dr0));
          theta_ph=atan(zph2/(r0+dr0));
          Tph = TSTARRAD*sqrt(RSTARRAD*4.65e-3/((r0+dr0)/cos(theta_ph))/5.2/2)
               /pow(RATIOPNUM/ABSFRACNUM,0.25);      // K
          csph = sqrt(Tph*GAMMA*k_b/mu/mp); // cm/s
          csph *= velocity_convert; // scale free
          Hph2 = csph/sqrt(G*MSTAR*MSTARRAD/pow(r0+dr0,3));  //scale free
          rhoph2 = sigma2/(2.5066*H2)*
                  exp(-zph2*zph2/2*(1/H2/H2-1/Hph2/Hph2)); 
          sigma_alt = sigma2/2*erf(zph2/1.414/H2) + rhoph2*sqrt(3.1415/2)*Hph2*(1-erf(zph2/1.414/Hph2));
          reduction = sigma2/2/(sigma_alt);
          if((r0+dr0)*theta>zph2){
            rhoph2*=reduction;
            rho2 = rhoph2*exp(-pow(((r0+dr0)*tan(theta))/(1.414*Hph2),2.0));
          } 
          else{
            sigma2*=reduction;
            rho2 = sigma2/(2.507*H2)*exp(-pow(((r0+dr0)*tan(theta))/(1.414*H2),2.0));
          }

          tau+=exp(0.5*log(rho1)+0.5*log(rho2));
          r0+=dr0;
        }
        tau*=opac*dr0/cos(theta);
      }
      else{
        //a simplified integration that doesn't result in Gammas
        //assuming Sigma=Sigma0 * (r/R0)^(-1) and simply H \propto r
        rho0 = (avg_rho_in/(2.507*H_azi[la]));
        if(S == 1.0) tau = opac * rho0*r*log(r/Si_subl_r); 
        else tau = opac*rho0*pow(r,S)/(1-S)*(pow(r,1-S)-pow(Si_subl_r,1-S));
        tau *= exp(-pow((r*tan(theta))/(1.414*H_azi[la]),2.0)) / cos(theta);            
      }    
      tau *= interior_factor;
      //check if the optical depth is too high/low or good enough
      if(tau < surf_depth){
        dtheta = 0.5*fabs(theta-theta_prev);
        theta_prev = theta;
        theta -= dtheta;
        trapped_count += 1;
        if(trapped_count>=200){
          masterprint("ERROR! inner edge not found i=%d\n",i);
          prs_exit(1);
        }
      }
      else if(tau > surf_depth+5.0e-4){
        dtheta = 0.5*fabs(theta-theta_prev);
        theta_prev = theta;
        theta += dtheta;
        trapped_count += 1;
      }
      else surf_found = 1;
    }

    ////////////////////////////////////////
    // step 1.5: setup theta & tau arrays //
    ////////////////////////////////////////

    //now that the surface is found, make arrays of theta vs tau
    //above the surface. Interpolate between these arrays when finding
    //interior optical depths for larger thetas. 
    surf_azi[la] = r*tan(theta)*5.2;
    theta_arr0[0] = theta;
    tau_arr0[0] = tau;
    for(ii=1;ii<nray;ii++){
      theta = theta_arr0[0]+(MAXTHETA-theta_arr0[0])*ii/(nray-1);
      if(use_gamma==1){
        //complete Sigma and H integration, resulting in Gamma functions
        c1 = opac*avg_rho_in*pow(r,S-F)/(2.5066*A);
        c2 = S+F+1;
        c3 = 2*F;
        c4 = tan(theta)*tan(theta)*pow(r,c3)/(2*A*A);        
        g1 = calc_upper_incomp_gamma((c2-1)/c3,c4*pow(r,-c3));
        g2 = calc_upper_incomp_gamma((c2-1)/c3,c4*pow(Si_subl_r,-c3));        
        tau = c1*pow(c4,(1-c2)/c3)/c3*(g1-g2);
      }
      else if(just_integrate==1){
        r0 = Si_subl_r;
        tau=0;
        //sigma1 = avg_rho_in*pow(r0/R0,avg_alpha_rho)*(1-r0/rpos[NGHY]) +
          //       dens_azi[la]*pow(r0/rpos[NGHY],avg_alpha_rho)*(r0/rpos[NGHY]);
        sigma1 = avg_rho_in*pow(r0/R0,avg_alpha_rho);
        H1 = ASPECTRATIO*pow(r0/R0,FLARINGINDEX)*r0*(1-r0/rpos[NGHY]) +
                 H_azi[la]*pow(r0/rpos[NGHY],FLARINGINDEX+1)*(r0/rpos[NGHY]);
        zph1 = get_photosphere(sigma1,H1,opacity->Rosseland,zph_depth,0.0)*sqrt(rpos[NGHY]/r0);
        theta_ph=atan(zph1/r0);
        Tph = TSTARRAD*sqrt(RSTARRAD*4.65e-3/(r0/cos(theta_ph))/5.2/2)
             /pow(RATIOPNUM/ABSFRACNUM,0.25);      // K
        csph = sqrt(Tph*GAMMA*k_b/mu/mp); // cm/s
        csph *= velocity_convert; // scale free
        Hph1 = csph/sqrt(G*MSTAR*MSTARRAD/pow(r0,3));  //scale free
        rhoph1 = sigma1/(2.5066*H1)*
                 exp(-zph1*zph1/2*(1/H1/H1-1/Hph1/Hph1)); 
        sigma_alt = sigma1/2*erf(zph1/1.414/H1) + rhoph1*sqrt(3.1415/2)*Hph1*(1-erf(zph1/1.414/Hph1));
        reduction = sigma1/2/(sigma_alt);
        if(r0*theta>zph1){
          rhoph1*=reduction;
          rho2 = rhoph1*exp(-pow((r0*tan(theta))/(1.414*Hph1),2.0));
        } 
        else{
          sigma1*=reduction;
          rho2 = sigma1/(2.507*H1)*exp(-pow((r0*tan(theta))/(1.414*H1),2.0));
        }
        for(jj=0;jj<nstep0;jj++){
          rho1 = rho2;
          //sigma2 = avg_rho_in*pow((r0+dr0)/R0,avg_rho_in)*(1-(r0+dr0)/rpos[NGHY]) +
            //       dens_azi[la]*pow((r0+dr0)/rpos[NGHY],avg_rho_in)*((r0+dr0)/rpos[NGHY]);
          sigma2 = avg_rho_in*pow((r0+dr0)/R0,avg_rho_in);
          H2 = ASPECTRATIO*pow((r0+dr0)/R0,FLARINGINDEX)*(r0+dr0)*(1-(r0+dr0)/rpos[NGHY]) +
                   H_azi[la]*pow((r0+dr0)/rpos[NGHY],FLARINGINDEX+1)*((r0+dr0)/rpos[NGHY]);
          zph2 = get_photosphere(sigma2,H2,opacity->Rosseland,zph_depth,0.0)* sqrt(rpos[NGHY]/(r0+dr0));
          theta_ph=atan(zph2/(r0+dr0));
          Tph = TSTARRAD*sqrt(RSTARRAD*4.65e-3/((r0+dr0)/cos(theta_ph))/5.2/2)
               /pow(RATIOPNUM/ABSFRACNUM,0.25);      // K
          csph = sqrt(Tph*GAMMA*k_b/mu/mp); // cm/s
          csph *= velocity_convert; // scale free
          Hph2 = csph/sqrt(G*MSTAR*MSTARRAD/pow((r0+dr0),3));  //scale free
          rhoph2 = sigma2/(2.5066*H2)*
                  exp(-zph2*zph2/2*(1/H2/H2-1/Hph2/Hph2)); 
          sigma_alt = sigma2/2*erf(zph2/1.414/H2) + rhoph2*sqrt(3.1415/2)*Hph2*(1-erf(zph2/1.414/Hph2));
          reduction = sigma2/2/(sigma_alt);
          if((r0+dr0)*theta>zph2){
            rhoph2*=reduction;
            rho2 = rhoph2*exp(-pow(((r0+dr0)*tan(theta))/(1.414*Hph2),2.0));
          } 
          else{
            sigma2*=reduction;
            rho2 = sigma2/(2.507*H2)*exp(-pow(((r0+dr0)*tan(theta))/(1.414*H2),2.0));
          }
          if(!isnan(rho1*rho2) && rho1*rho2>0.0) tau+=exp(0.5*log(rho1)+0.5*log(rho2));
          r0+=dr0;
        }
        tau*=opac*dr0/cos(theta);
      }
      else{
        //a simplified integration that doesn't result in Gammas
        //assuming Sigma=Sigma0 * (r/R0)^(-1) and simply H \propto r
        rho0 = (avg_rho_in/(2.507*H_azi[la]));
        if(S == 1.0) tau = opac * rho0*r*log(r/Si_subl_r); 
        else tau=opac*rho0*pow(r,S)/(1-S)*(pow(r,1-S)-pow(Si_subl_r,1-S));
        tau *=exp(-pow((r*tan(theta))/(1.414*H_azi[la]),2.0)) / cos(theta);  
      }    
      tau *= interior_factor;       
      theta_arr0[ii] = theta;  //arrays of theta and tau that act as
      tau_arr0[ii] = tau;      //the starting points for the integration
    }
  

  for(i=0;i<size_x;i++){
    //////////////////////////////////////////////
    //// step 2: find the rest of the surface ////
    //////////////////////////////////////////////
    for(ii=0;ii<nray;ii++){
      tau_arr[ii] = tau_arr0[ii];
      theta_arr[ii] = theta_arr0[ii];
    }

    ii=0; //index of the lowest angle ray that still needs to be calculated
    j=j_start-1;
    sigma2 = dens_azi[la];
    H2 = H_azi[la];
    zph2 = get_photosphere(sigma2,H2,opacity->Rosseland,zph_depth,0.0);
    theta_ph=atan(zph2/rpos[j]);
    Tph = TSTARRAD*sqrt(RSTARRAD*4.65e-3/(rpos[j]/cos(theta_ph))/5.2/2)
               /pow(RATIOPNUM/ABSFRACNUM,0.25);      // K
    csph = sqrt(Tph*GAMMA*k_b/mu/mp); // cm/s
    csph *= velocity_convert; // scale free
    Hph2 = csph/sqrt(G*MSTAR*MSTARRAD/pow(rpos[j]*rpos[j]+zph2*zph2,1.5));  //scale free
    rhoph2 = sigma2/(2.5066*H2)*
                  exp(-zph2*zph2/2*(1/H2/H2-1/Hph2/Hph2)); 
    sigma_alt = sigma2/2*erf(zph2/1.414/H2) + rhoph2*sqrt(3.1415/2)*Hph2*(1-erf(zph2/1.414/Hph2));
    reduction = sigma2/2/(sigma_alt);
    rhoph2 *= reduction;
    sigma2*=reduction;
    for(j=j_start;j<j_end;j++){
      //first, calculate the optical depth at this new radial index for all
      //of the rays that still need to be calculated
      r = rpos[j];
      dr = rpos[j]-rpos[j-1];
      
      H1 = H2;
      zph1 = zph2;   
      Hph1 = Hph2;  //scale free
      rhoph1 = rhoph2;
      sigma1 = sigma2;

      sigma2=dens_azi[la];
      H2 = H_azi[la];
      zph2 = get_photosphere(sigma2,H2,opacity->Rosseland,zph_depth,zph1);
      if(theta_ph<atan(zph2/rpos[j])) theta_ph=atan(zph2/rpos[j]) ; 
      Tph = TSTARRAD*sqrt(RSTARRAD*4.65e-3/(rpos[j]/cos(theta_ph))/5.2/2)
               /pow(RATIOPNUM/ABSFRACNUM,0.25);      // K
      csph = sqrt(Tph*GAMMA*k_b/mu/mp); // cm/s
      csph *= velocity_convert; // scale free
      Hph2 = csph/sqrt(G*MSTAR*MSTARRAD/pow(rpos[j]*rpos[j]+zph2*zph2,1.5));  //scale free
      rhoph2 = sigma2/(2.5066*H2)*
                  exp(-zph2*zph2/2*(1/H2/H2-1/Hph2/Hph2)); 
      sigma_alt = sigma2/2*erf(zph2/1.414/H2) + rhoph2*sqrt(3.1415/2)*Hph2*(1-erf(zph2/1.414/Hph2));
      reduction = sigma2/2/(sigma_alt);
      rhoph2 *= reduction;
      sigma2*=reduction;

      for(i2=ii;i2<nray;i2++){  
        theta = theta_arr[i2];
        dl = dr/cos(theta);
        z = (r-dr)*tan(theta);
        //tau is calc'd with log interpolation of density between points
        if(z<zph1){
          //solving for log(density) since we use a log interpolation
          rho1 =  log((sigma1/(2.507*H1)))
                 -pow(z/(1.414*H1),2.0); 
        }
        else rho1 = log(rhoph1)-pow(z/(1.414*Hph1),2.0);
        z = r*tan(theta);
        if(z<zph2){
           //solving for log(density) since we use a log interpolation
          rho2 = log(sigma2/(2.507*H2))
                  -pow(z/(1.414*H2),2.0);
        }
        else rho2 = log(rhoph2)-pow(z/(1.414*Hph2),2.0);

        tau_arr[i2]+=opac*exp((rho1+rho2)/2.0)*dl;
      }
      while(tau_arr[ii+1] > surf_depth){
        ii++; //once the next lowest ray is >2/3, we don't need the lowest
      }   
      //use a gausian-like interpolation to find the tangent
      //of the angle at which tau=2/3 for this radial index
      c1 = (log(surf_depth)-log(tau_arr[ii]))
          /(log(tau_arr[ii+1])-log(tau_arr[ii]));
      c2 = c1*(tan(theta_arr[ii+1])*tan(theta_arr[ii+1])
           - tan(theta_arr[ii])*tan(theta_arr[ii]));
      tan_theta = sqrt(c2+tan(theta_arr[ii])*tan(theta_arr[ii]));
      //record the surface height
      surf_azi[la] = r*tan_theta*5.2; 
    } 
    //adjust the maximum theta value based on the current back edge
    if(atan(tan_theta)>new_max_theta) new_max_theta=atan(tan_theta);   
  }
  //set the new maximum theta value to be used for the next timestep
  MAXTHETA = 1.05*new_max_theta; 
  
  //extend the curvature of the disk to get surface of inner ghost cells(+1)
  double dz2,dz1,ddz; 
  for(i=0;i<size_x;i++){  
    dz2 = surf_azi[size_x*(NGHY+3) + i]-surf_azi[size_x*(NGHY+2) + i];
    dz1 = surf_azi[size_x*(NGHY+2) + i]-surf_azi[size_x*(NGHY+1) + i];
    ddz = dz2-dz1;
    for(j=NGHY;j>-1;j--){
      surf_azi[la] =  surf_azi[size_x*(j+1) + i]-dz1-ddz*(j-NGHY+1);
    }
    dz2 = surf_azi[size_x*(size_y-NGHY-2) + i]-surf_azi[size_x*(size_y-NGHY-3) + i];
    dz1 = surf_azi[size_x*(size_y-NGHY-1) + i]-surf_azi[size_x*(size_y-NGHY-2) + i];
    ddz = dz1-dz2;
    for(j=size_y-NGHY;j<size_y;j++){
      surf_azi[la] =  surf_azi[size_x*(j-1) + i]+dz1+ddz*(j-(size_y-NGHY-1));
    }
  }

  free(tau_arr);
  free(theta_arr);
  free(tau_arr0);
  free(theta_arr0);

}

double get_photosphere(double sigma, double H, double opac, double depth, float limit){
  double z, tau_v, tau, dz, rho0;
  double tolerance, maclaurin, x,tt1, tt2;
  double H_root2 = H*1.414;
  int photosphere_found=0;
  int count = 0;
  tolerance = depth*2.0e-2;
  tau_v = sigma*opac/2.0;
  z = 0.0;
  if(tau_v>depth){
    x = (tau_v - depth)*2/opac/sigma;
    x = 1-x*x;        
    tt1 = 2/(M_PI*0.147) + 0.5 * log(x);
    tt2 = 1/(0.147) * log(x);
    z = H_root2*sqrt(-tt1 + sqrt(tt1*tt1 - tt2));
  }
  if(z<limit) z=limit;

  return z;
}


double calc_upper_incomp_gamma(double s, double x){
  //a function to calculate the upper incomplete gamma function
  //by subtracting the lower incomplete gamma from the full gamma
  double gamma, gamma0, gamma_i, num, denom, val;
  int k;
  denom = s;
  gamma0 = pow(x,s)*exp(-x)/s;
  gamma = gamma0;
  gamma_i = gamma0;
   
  k=1;
  while(gamma_i > gamma0*1e-4){
    denom *= (s+k);
    num = pow(x,s)*exp(-x)*pow(x,k);
    gamma_i = num/denom;
    gamma += gamma_i;
    k++;
  }
  
  val = tgamma(s)-gamma;
  
  return val;
}


void calc_surf_lookup(struct disk_parameters *DP, struct disk_opacity *opacity,
                   double *dens_azi, double *H_azi, double *surf_azi,int opac_i,
                   int size_x,int size_y,double avg_rho_in, double *rpos)
{
  //A function to calculate the surface of the disk. The optical depth is 
  //calculated along fixed rays and the exact height of the tau=2/3 surface 
  //for the starlight is found by interpolation at the radial grid points
  double rho0, rho1, rho2, tau, dtheta;
  double theta, theta_prev, tan_theta; 
  double opac, dl, Si_subl_r;
  int ii,i,i2,j,k,j2;
  int j_shift=(Stride*CPU_Rank/Nx-2*NGHY*CPU_Rank);
  int surf_found, trapped_count,j_start,j_end,n;
  double new_max_theta=0.0;
  double A,S,F, c1, c2, c3, c4, g1, g2, t_prev, s_prev;
  double dr;  //scale-free
  double r0;  //scale-free (+.5 bc Ymed was used)
  double zph1, zph2, csph, rhoph1, rhoph2, rho2o, Hph1, Hph2, Tph, z, theta_ph;
  double sigma_alt, reduction;
  double energy_convert = G/G_CGS*pow(MSTAR/MSTAR_CGS,2)*pow(R0_CGS/R0,3);
  double velocity_convert = pow(R0_CGS/R0*G/G_CGS*MSTAR/MSTAR_CGS,0.5);
  real mu = 2.3;
  real mp = 1.67e-24 ;  //g
  real k_b = 1.381e-16; // erg K-1
  real gamma = 1.667;
  real StfBolt = 5.67e-5; //erg cm-2 s-2 K-4
  double abs_frac = opacity->abs_frac;
  
  //parameters used to extrapolate the density inside the grid inner edge
  S = SIGMASLOPE;
  F = FLARINGINDEX;
  //rough estimate of Silica sublimation radius
  Si_subl_r = DP->T_star*DP->T_star/1250.0/1250.0*DP->R_star*R0/R0_CGS;
  if(YMIN<Si_subl_r){
    printf("Error: inner edge (%f) < estimate for Si sublimation (%f)\n",
           YMIN,Si_subl_r);
    prs_exit(EXIT_FAILURE);
  }
  double dr0;
  double sigma1, sigma2, H1, H2;
  int just_integrate=0;
  int use_gamma=0;
  int jj;
  int nstep0=50;
  dr0 = (rpos[NGHY]-Si_subl_r)/nstep0;


  real omega, r;

  int nray= 120;//y_tot/12;
  double theta_factor;
  double *theta_arr;
  double *tau_arr;
  theta_arr = (real *) malloc(sizeof(double *)*nray);
  tau_arr = (real *) malloc(sizeof(double *)*nray);


  double interior_factor=0.02;

  //opacity already scale free
  opac = (opacity->Rosseland);

  FILE *fsurf = fopen("/Users/leiendeck/research/fargo_grainsize/surface_hd163296.txt","r");
  if(fsurf==NULL) printf("error opening file  ");
  double surf_arr[256];
  for (i=0; i < 256; i++) {
    if(fscanf(fsurf,"%lf",&surf_arr[i])!=1) printf("error!  ");
  }
  fclose(fsurf);
  double dh = 0.20059899/5.2, dhrem;
  int kz;


  for(i=0;i<size_x;i++){
    for(j=NGHY;j<NY+NGHY;j++){
      surf_azi[la] = surf_arr[j-NGHY];
    }
  }
  double dz2,dz1,ddz; 
  for(i=0;i<size_x;i++){  
    dz2 = surf_azi[size_x*(NGHY+3) + i]-surf_azi[size_x*(NGHY+2) + i];
    dz1 = surf_azi[size_x*(NGHY+2) + i]-surf_azi[size_x*(NGHY+1) + i];
    ddz = dz2-dz1;
    for(j=NGHY;j>-1;j--){
      surf_azi[la] =  surf_azi[size_x*(j+1) + i]-dz1-ddz*(j-NGHY+1);
    }
    dz2 = surf_azi[size_x*(size_y-NGHY-2) + i]-surf_azi[size_x*(size_y-NGHY-3) + i];
    dz1 = surf_azi[size_x*(size_y-NGHY-1) + i]-surf_azi[size_x*(size_y-NGHY-2) + i];
    ddz = dz1-dz2;
    for(j=size_y-NGHY;j<size_y;j++){
      surf_azi[la] =  surf_azi[size_x*(j-1) + i]+dz1+ddz*(j-(size_y-NGHY-1));
    }
  }


  free(tau_arr);
  free(theta_arr);

 
  
}


void output_density_struct(struct disk_parameters *DP, struct disk_opacity *opacity,
  double *dens_azi, double *H_azi, double *surf_azi,int opac_i,
  int size_x,int size_y, double *rpos)
{
//A function to calculate the surface of the disk. The optical depth is 
//calculated along fixed rays and the exact height of the tau=2/3 surface 
//for the starlight is found by interpolation at the radial grid points
double rho0, rho1, rho2, tau, dtheta,r,rsphere;
double theta, theta_prev, tan_theta; 
int ii,i,i2,j,k,j2;
double dr, r0, c, rcyl;  
double zph, zph_real, csph, rhoph, H, Hph, Tph, z, theta_ph, theta_max;
double rho, sigma;
double sigma_alt, reduction;
double energy_convert = G/G_CGS*pow(MSTAR/MSTAR_CGS,2)*pow(R0_CGS/R0,3);
double velocity_convert = pow(R0_CGS/R0*G/G_CGS*MSTAR/MSTAR_CGS,0.5);
double density_convert = MSTAR/MSTAR_CGS/pow(R0/R0_CGS,3);
real mu = 2.3;
real mp = 1.67e-24 ;  //g
real k_b = 1.381e-16; // erg K-1
real gamma = 1.667;
real StfBolt = 5.67e-5; //erg cm-2 s-2 K-4
double zph_depth = 0.1; //1.0/3; //3e-4 

if(CPU_Number>1){
printf("Trying to output density structure with >1 CPU\n");
exit(0);
}

printf("going to print the 3D density structure somewhere...\n");

FILE *fp3;
char str3[180];
sprintf(str3,"/Users/leiendeck/research/radmc3d-2.0-master/comparison_tests/spiral_test/production/fargo_shadow/3D_dens_cyl.txt");
fp3 = fopen(str3, "w+");
FILE *fp2;
char str2[180];
sprintf(str2,"/Users/leiendeck/research/radmc3d-2.0-master/comparison_tests/spiral_test/production/fargo_shadow/photosphere.txt");
fp2 = fopen(str2, "w+");

r0 = rpos[0];
theta_max = 0.32;
for(i=0;i<size_x;i++){
for(ii=0; ii<48; ii++){
//for(k=0;k<128;k++){
//z = k*3.5/5.2/128;
theta = theta_max*ii/48;
zph = 0.0;
for(j=NGHY;j<size_y-NGHY;j++){
rsphere = rpos[j];
rcyl = rsphere*cos(theta);
z = rsphere*sin(theta);   
//theta = atan(z/rcyl);     

/*if(((toupper(*SPACING)) == 'L') && ((toupper(*(SPACING+1))) == 'O')){
dr = (log(rpos[NY+2*NGHY-1]*5.2)-log(rpos[0]*5.2))/(NY+2*NGHY-1);
j2 = (int)((log(rsphere)-log(r0))/dr);
if(j2<0) printf("r=%f YMIN=%f dr=%E\n",rsphere,YMIN,dr);
if(j2<size_y-1) dr = (rpos[j2+1]-rpos[j2]);
else  dr = (rpos[j2]-rpos[j2-1]);
}*/
//else{
dr = (rpos[1]-rpos[0]);
j2 = (int) ((rcyl-r0)/dr);
//} 
if(j2<0) j2 = 0;
c = (rcyl-rpos[j2])/dr;
if(c<0) c=0.0;

rcyl = (1.0-c)*(rpos[j2]) + c*(rpos[(j2+1)]);
sigma = (1.0-c)*(dens_azi[i+j2*size_x]) + c*(dens_azi[i+(j2+1)*size_x]);
H = (1.0-c)*(H_azi[i+j2*size_x]) + c*(H_azi[i+(j2+1)*size_x]);
zph = get_photosphere(sigma,H,opacity->Rosseland,zph_depth,zph);
zph_real = get_photosphere(sigma,H,opacity->Rosseland,2.0/3,0.0);
if(j==NGHY) theta_ph=atan(zph/rcyl);
else if(theta_ph<atan(zph/rcyl)) theta_ph=atan(zph/rcyl) ; 
Tph = TSTARRAD*sqrt(RSTARRAD*4.65e-3/(rcyl/cos(theta_ph))/5.2/2)
/pow(RATIOPNUM/ABSFRACNUM,0.25);      // K
csph = sqrt(Tph*GAMMA*k_b/mu/mp); // cm/s
csph *= velocity_convert; // scale free
Hph = csph/sqrt(G*MSTAR*MSTARRAD/pow(rcyl*rcyl+zph*zph,1.5));  //scale free
rhoph = sigma/(2.5066*H)*
   exp(-zph*zph/2*(1/H/H-1/Hph/Hph)); 
sigma_alt = sigma/2*erf(zph/1.414/H) + rhoph*sqrt(3.1415/2)*Hph*(1-erf(zph/1.414/Hph));

reduction = sigma/2/(sigma_alt);
rhoph *= reduction;
sigma *= reduction;
if(z<zph) rho =  (sigma/(2.507*H))*exp(-pow(z/(1.414*H),2.0)); 
else rho = rhoph*exp(-pow(z/(1.414*Hph),2.0));
rho *= 1/density_convert;

fprintf(fp3, "%.4E \n",rho);
if(ii==0) fprintf(fp2, "%.4E \n",zph_real*5.2);
if((i==0 && j==67) || (i==0 && ii==0)) printf("j=%d theta=%.4f rsph=%.4f rcyl=%.4f  z=%f j2=%d c=%.3f   -   rho=%.3E   zph=%f vs %f \n",
j,theta,rsphere,rcyl,z,j2,c,rho*density_convert,zph,zph_real);
}
//fprintf(fp3,"\n");
}   
}
fclose(fp3);
fclose(fp2);
}
