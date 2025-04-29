#include "fargo3d.h"

void _CondInit(real epsilon) {
  
  int i,j,k;
  real r, omega;
  
  real *rho  = Density->field_cpu;
  real soundspeed;
#ifdef ADIABATIC
  real *e   = Energy->field_cpu;
#endif
#ifdef ISOTHERMAL
  real *cs   = Energy->field_cpu;
#endif
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
#ifdef STOKES2POP
  real *grainsize = GrainSize->field_cpu; //HSL
#endif
#ifdef RUNRADTRANS
  real *H =HTherm->field_cpu; //HSL
#endif
  
  real rhog, rhod;
  real vk, H2;
  
  i = j = k = 0;
  
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
	
	      r     = Ymed(j);
	      omega = sqrt(G*MSTAR/r/r/r);                       //Keplerian frequency
        //gas surface density
	      rhog  = SIGMA0*pow(r/R0,-SIGMASLOPE);    
        rhod  = rhog*epsilon;                           

        if (Fluidtype == GAS) {
          rho[l]   = rhog;
          vphi[l]  = omega*r*sqrt(1.0 + pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
                (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
          vr[l]    = 0.0 ;
          soundspeed  = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;
#ifdef ISOTHERMAL
          cs[l] = soundspeed;
#endif
#ifdef ADIABATIC
          e[l] = pow(soundspeed,2)*rho[l]/(GAMMA-1.0);
#endif
#ifdef STOKES2POP          
          grainsize[l] = 0.0; //HSL
#endif
#ifdef RUNRADTRANS
          H2 = soundspeed/omega; 
          H[l] = H2;

#endif
	}
	
	if (Fluidtype == DUST) {
	  rho[l]  = rhod;
	  vphi[l] = omega*r;
	  vr[l]   = 0.0; 
#ifdef ISOTHERMAL
          cs[l] = 0.0;
#endif
#ifdef ADIABATIC
          e[l] = 0.0; 
#endif
#ifdef STOKES2POP          
          grainsize[l] = 2.0; //HSL
#endif
#ifdef RUNRADTRANS
          H[l] = H2;
#endif
	}
	
	vphi[l] -= OMEGAFRAME*r;
	
      }
    }
  }
}

void CondInit() {
  int id_gas = 0;
  int feedback = YES;
  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit(EPSILON1);

  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;
#ifdef STOKES2POP
  int i,j,k; //added by harrison
  real *asize; //HSL
  //real *alpha_2pop;
  real dust_size; //HSL
  //real *e_grow; //HSL
  real rho_s = RHO_S*pow(R0_CGS/R0,3)*MSTAR/MSTAR_CGS;

  real *rho;  
  rho  = Fluids[0]->Density->field_cpu;
#endif
  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
#ifdef STOKES2POP
    if (id_dust==1) _CondInit(EPSILON1);
    if (id_dust==2) _CondInit(EPSILON2);

    /*all added by harrison*/
    i = j = k = 0;
    asize = Fluids[id_dust]->GrainSize->field_cpu;
    //alpha_2pop=Fluids[id_dust]->Alpha0->field_cpu;
    //e_grow = Fluids[id_dust]->E_grow->field_cpu;

    if (id_dust==1)  dust_size = DUST_A1;
    if (id_dust==2)  dust_size = DUST_A2;
    dust_size*=R0/R0_CGS;

    for (k=0; k<Nz+2*NGHZ; k++) {
      for (j=0; j<Ny+2*NGHY; j++) {
        //printf("check initial alpha:id=%d j=%d %E \n",id_dust,j,rho[l]);
        for (i=0; i<Nx+2*NGHX; i++) {
          //original
          //set the grain sizes
          if (id_dust==1) asize[l] = DUST_A1*R0/R0_CGS;
          if (id_dust==2) asize[l] = DUST_A2*R0/R0_CGS;
        }
      }
    }
#else
    _CondInit(EPSILON1);
#endif

  }

#ifndef STOKES2POP
  ColRate(INVSTOKES1, id_gas, 1, feedback);
  ColRate(INVSTOKES2, id_gas, 2, feedback);
#endif

}

