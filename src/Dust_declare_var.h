real r;
real r_p1;
real r_m1;
real P;
real P_p1;
real P_m1;
real omega_p1;
real omega_m1;
real gamma;
int n;

real alpha_rn;  //added by HSL
real a0;        // grain size small pop
real a1;        //grain size large pop
real a1i;        //grain size large pop
real rho_s = RHO_S*pow(R0_CGS/R0,3)*MSTAR/MSTAR_CGS;//2.367e8; //dust material density 1.5 gm cm-3 
                   // x (5.2 au in cm)^3/M_sun (g) to become scale free
real mu = 2.3;
real m_p = 1.67e-24*MSTAR/MSTAR_CGS;  //proton mass g->scale free

real a_fr_ep;
real a_fr_st;
real a_fr;
real a_dr;
real a_df;
real a_St1;
real a_grow;

real InvSt1, InvSt2;

real fudge_fr = 0.37;  //fudge factor
real fudge_dr = 0.55; 
real v_frag = VFRAG*sqrt(G/G_CGS*MSTAR/MSTAR_CGS/R0*R0_CGS); //cm s-1->scl free
real sig_h2 = 2e-15*R0*R0/R0_CGS/R0_CGS;  //collsion cross sec cm^2->scl free
real n_d, t_col;
real v_turb11;
real dvy;

real tau_grow;
real tau_grow_ph;
real E_stick = 1;

real mass_sweep;
real mass_frag;
real v_sign;
real dv_sweep;
real dv_frag;
real v11;
real v01;
real sig01;
real sig11;
real H_0; //HSL
real H_1; 
real F;
real plaw;
//real vg, Re, St1, St2, tf2, tf1, t0, ts, beta, St_rat;


//declerations for ormel & cuzzi 2007 turbulent relative velocity
real c0 = 1.6015125;
real c1 = -0.63119577;
real c2 = 0.32938936;
real c3 = -0.29847604;
real ya = 1.6;
real yap1inv = 1.0 / (1.0 + ya);
real OmKinv, Re, ReInvSqrt, vn, vs, ts, vg2;
real St1, St2, StL, StS, StM, epsLM, epsLS, tauL, tauS,tauM , ys, h1LS, h1LM , h2LS, h2LM;


//real a_grow_1[size_y*size_x];
//real a_grow_2[size_y*size_x];
