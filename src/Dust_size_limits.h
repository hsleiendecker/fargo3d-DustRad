/* 

new file by HSL
trying to update the alpha values (inv stokes number) based on 2pop

*/

// In the implementation below, alpha --> 1/St
// using NFLUID = 1 as the small population
// and NFLUID = 2 as the large population

InvSt1 = 2*rho[0][l]/(DUST_A1/R0_CGS*R0*rho_s*3.1415);
InvSt2 = 2*rho[0][l]/(asize[2][l]*rho_s*3.1415);

if(id2==1){

  a0 = DUST_A1/R0_CGS*R0;
  a1 = asize[2][l];
  a1i=a1;

  r = (id1+id3)*ymed(j)+id2*ymin(j);
  omega = sqrt(G*MSTAR/(r*r*r));

#ifdef RADTRANSFER
  H_0 = H_gas[l];
#else
  H_0 = cs[l]/omega;
#endif
  if(1.0/InvSt2 < 0.5){
    H_1 = H_0 * pow(ALPHA/(1.0/InvSt2
              *(1+1.0/(InvSt2*InvSt2))),0.5);
  }
  else{
    H_1 = H_0*pow(ALPHA/(2.0*(1+1.0/(InvSt2*InvSt2))),0.5);
  }
  if(H_1 > H_0) H_1 = H_0;

    P = rho[0][l]*omega * cs[l];
  if(j > 0){
    r_m1 = (id1+id3)*ymed(j-1)+id2*ymin(j-1);
    omega_m1 = sqrt(G*MSTAR/(r_m1*r_m1*r_m1));
    P_m1 = rho[0][lym]*omega_m1 * cs[lym];
  }
  if(j < size_y-1){
    r_p1 = (id1+id3)*ymed(j+1)+id2*ymin(j+1);
    omega_p1 = sqrt(G*MSTAR/(r_p1*r_p1*r_p1));
      P_p1 = rho[0][lyp]*omega_p1 * cs[lyp];
  }

  if(j < 1){
    gamma = r/P*(P_p1-P)/(r_p1-r);
  }
  else if(j == size_y-1){
    gamma = r/P*(P-P_m1)/(r-r_m1);
  }
  else{
    gamma = r/P*(P_p1-P_m1)/(r_p1-r_m1);
  }

  // start with a_fr since that has no gradient dependence
  a_fr_ep = fudge_fr*rho[0][l]*v_frag*v_frag
             /(4.712*ALPHA*rho_s * cs[l] * cs[l]);

  a_fr = a_fr_ep;

  if (1/InvSt2 > 0.5){
    a_fr_st = 0.691*pow(0.5/(sig_h2/(2.5066 * cs[l]/omega*mu*m_p))/(ALPHA*rho_s),0.5)
              *v_frag/cs[l];
    if(a_fr_st < a_fr){
      a_fr = a_fr_st;
    }
  }

  double eta;

  eta = GAMMA*P/r/(2*rho[0][l]/(2.506*cs[l]/omega)*r*omega*omega);

  //E_drift=1 from two-pop-py
  //version from two-pop-py
  /*a_dr = E_stick*fudge_dr*0.6366*(rho[1][l]+rho[2][l])/rho_s*r*r*omega*omega
         /(fabs(gamma) * cs[l] * cs[l]);
    a_df = fudge_fr*2.0*rho[0][l]/(3.142*rho_s)*v_frag*pow(G*MSTAR/r,0.5)
         /(fabs(gamma) * cs[l] * cs[l]*(1-0.5)); */
  //version from drazkowska et al 2019
  a_dr = fudge_dr*(rho[2][l]/0.97)/(2*fabs(eta)) * 2/(3.1415*rho_s);
  //Drazkowska et al 2019 version
  a_df = fudge_fr*v_frag/(fabs(eta)*(r*omega)) * 0.6366*rho[0][l]/rho_s;

  if(a0 > a_df) a_df = a0;

  /*assume erosion limits particles to St=1 and treat like frag*/
  a_St1 = 0.6366*rho[0][l]/rho_s; 
  if(a_St1 < a_fr) a_fr = a_St1;
  if(a_St1 < a_dr) a_dr = a_St1;

  //calc the growth timescale and thus a_1(t)
  //tau_grow_ph appears to have units of mass/area/time
  tau_grow_ph = 1e-100*pow(R0_CGS/R0,2)*MSTAR/MSTAR_CGS/SCALE_FREE_T;
  if(E_stick*(rho[1][l]+rho[2][l])*omega > tau_grow_ph){
    tau_grow_ph = E_stick*(rho[1][l]+rho[2][l])*omega ;
  }
  tau_grow = rho[0][l]/tau_grow_ph;
        
  tau_grow_ph = dt/tau_grow;  
  a_grow = a1*exp(tau_grow_ph);
   
  a1 = a_fr;
  plaw = 0.25/0.75;
  if(a_dr < a1){
    a1 = a_dr;
    plaw = 0.03/.97;
  }
  if(a_df < a1){
    a1 = a_df;
  }
  if(a0 > a1){  
    a1 = a0;    
  }
  if(a_grow < a1){
    a1 = a_grow;
  }

  /*if(a1 < a1i && PhysicalTime>2.0){
    if (a1<a1i*exp(-20*tau_grow_ph)){
      a1 = a1i*exp(-20*tau_grow_ph);
    }
  }*/

  if(a1<a0) a1 = a0;

  asize[2][l]=  a1;
  asize[1][l]= DUST_A1/R0_CGS*R0;
}
