//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void _collisions_cpu(real dt, int id1, int id2, int id3, int option) {
  
//<USER_DEFINED>

  real *rho[NFLUIDS];
  real *velocities_input[NFLUIDS];
  real *velocities_output[NFLUIDS];
#ifdef STOKES2POP
  real *asize[NFLUIDS]; //HSL
#endif
#ifdef RADTRANSFER
  real *H_gas; //HSL
#endif

  int ii;
#ifdef RADTRANSFER
  H_gas = Fluids[0]->HTherm->field_CPU;//HSL
#endif

  for (ii=0; ii<NFLUIDS; ii++) {

    INPUT(Fluids[ii]->Density);
    rho[ii]  = Fluids[ii]->Density->field_cpu;
#ifdef STOKES2POP
    asize[ii] = Fluids[ii]->GrainSize->field_cpu; //HSL
#endif

    //Collisions along X
    #ifdef X
    if (id1 == 1) {
      if (option == 1) {
	INPUT(Fluids[ii]->Vx_temp);
	OUTPUT(Fluids[ii]->Vx_temp);
	velocities_input[ii] = Fluids[ii]->Vx_temp->field_cpu;
	velocities_output[ii] = Fluids[ii]->Vx_temp->field_cpu;
      }
      if (option == 0) {
	INPUT(Fluids[ii]->Vx);
	OUTPUT(Fluids[ii]->Vx_half);
	velocities_input[ii] = Fluids[ii]->Vx->field_cpu;
	velocities_output[ii] = Fluids[ii]->Vx_half->field_cpu;
      }
    }
    #endif
    
    //Collisions along Y
    #ifdef Y
    if (id2 == 1) {
      if (option == 1) {
	INPUT(Fluids[ii]->Vy_temp);
	OUTPUT(Fluids[ii]->Vy_temp);
	velocities_input[ii] = Fluids[ii]->Vy_temp->field_cpu;
	velocities_output[ii] = Fluids[ii]->Vy_temp->field_cpu;
      }
      if (option == 0) {
	INPUT(Fluids[ii]->Vy);
	OUTPUT(Fluids[ii]->Vy_half);
	velocities_input[ii] = Fluids[ii]->Vy->field_cpu;
	velocities_output[ii] = Fluids[ii]->Vy_half->field_cpu;
      }
    }
    #endif
    
    //Collisions along Z
    #ifdef Z
    if (id3 == 1) {
      if (option == 1) {
	INPUT(Fluids[ii]->Vz_temp);
	OUTPUT(Fluids[ii]->Vz_temp);
	velocities_input[ii] = Fluids[ii]->Vz_temp->field_cpu;
	velocities_output[ii] = Fluids[ii]->Vz_temp->field_cpu;
      }
      if (option == 0) {
	INPUT(Fluids[ii]->Vz);
	OUTPUT(Fluids[ii]->Vz_half);
	velocities_input[ii] = Fluids[ii]->Vz->field_cpu;
	velocities_output[ii] = Fluids[ii]->Vz_half->field_cpu;
      }
    }
    #endif
  }

// External
  int size_x = XIP; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
// Internal
  int i;
  int j;
  int k;
#ifdef STOKES2POP
#ifdef ISOTHERMAL
  real *cs;         //HSL
  cs = Fluids[0]->Energy->field_cpu; //HSL
#endif
#ifdef ADIABATIC
  real *cs;
  cs = (real*)malloc(size_x*size_y*size_z*sizeof(real));
  real *energy_dens;
  energy_dens = Fluids[0]->Energy->field_cpu;
  i = j = k = 0;
#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=XIM; i<size_x; i++) {
#endif
        cs[l] = pow((GAMMA-1.0)*energy_dens[l]/rho[0][l],0.5);
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif

#endif
#endif



//<EXTERNAL>
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real* alpha = Alpha;
//<\EXTERNAL>

//<INTERNAL>
  int ll;
  int o;
  int p;
  int q;
  int ir;
  int ir2;
  int ir_max;
  int ic;
  real max_value;
  real factor;
  real big;
  real temp;
  real sum;
  int idm;
  real b[NFLUIDS];
  real m[NFLUIDS*NFLUIDS];  
  real omega;
  real rho_p;
  real rho_o;
  real rho_q; 

#ifdef STOKES2POP 
#include "Dust_declare_var.h"
#endif
//<\INTERNAL>

//<CONSTANT>
// real Alpha(NFLUIDS*NFLUIDS);
//<\CONSTANT>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=XIM; i<size_x; i++) {
#endif
//<#>

#ifdef STOKES2POP
#include  "Dust_size_limits.h" //HSL
#ifdef DUSTMASSTRANSFER
#include "Dust_mass_transfer.h"
#endif
#endif
	
#include  "collision_kernel.h"
#include  "gauss.h"
	
#ifdef COUPLEDUSTS
        /*forces the large and small dust population to drift together using
          their mass weighted velocity.  - HSL*/
        velocities_output[0][l] = b[0];
        velocities_output[1][l] = (b[1]*rho[1][l]+b[0]*rho[0][l])/(rho[1][l]+rho[0][l]);
        velocities_output[2][l] = (b[1]*rho[1][l]+b[2]*rho[2][l])/(rho[1][l]+rho[2][l]);
#else
	      for (o=0; o<NFLUIDS; o++) {
	        velocities_output[o][l] = b[o];
	      }    
#endif

//<\#>
#ifdef X
      }
      /*if(id2==1){
      fprintf(fp_f,"\n");
      fprintf(fp_d,"\n");
      fprintf(fp_df,"\n");
      fprintf(fp_g,"\n");
      }*/
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
  
  /*fclose(fp_f);
  fclose(fp_d);
  fclose(fp_df);
  fclose(fp_g);
  if(id2==1){
  for(i=0;i<1e10;i++){
    printf("*");
  }
  }*/

#ifdef STOKES2POP
#ifdef ADIABATIC
  free(cs);
#endif
#endif
//<\MAIN_LOOP>
}

void Collisions(real dt, int option) {

  //Input and output velocities are the same Fields
  if (option == 1) {
#ifdef X
    //Collisions along the X direction
    FARGO_SAFE(_collisions(dt,1,0,0,option));
#endif
#ifdef Y
    //Collisions along the Y direction
    FARGO_SAFE(_collisions(dt,0,1,0,option));
#endif
#ifdef Z
    //Collisions along the Z direction
    FARGO_SAFE(_collisions(dt,0,0,1,option));
#endif
  }
  
  //Input and output velocities are not the same Fields
  if (option == 0) {
#ifdef X
    //Collisions along the X direction
    FARGO_SAFE(_collisions(dt,1,0,0,option));
#endif
#ifdef Y
    //Collisions along the Y direction
    FARGO_SAFE(_collisions(dt,0,1,0,option));
#endif
#ifdef Z
    //Collisions along the Z direction
    FARGO_SAFE(_collisions(dt,0,0,1,option));
#endif
  }
}
