
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <mpi.h>

#define TAU_THICK (2.0/3.0)

ex void RadTransfer(real dt); //HSL
ex void (_radtransfer_azi_cpu)(real dt); //HSL

ex void DustEnergyCorrection();

ex void thermal_relaxation(real dt);


/* fundamental and astronomical constants */
extern const double 
  stefan_boltzmann,
  k_boltzmann,
  meanmolwt,
  gamma_th,
  AU_in_cm,
  year,
  grav_const,
  solar_mass,
  earth_mass;

/* structure that contains the disk/planet model parameters */
struct disk_parameters {
  double M_star, R_star, T_star, M_dot;
  double a_au, M_planet;
  double alpha, mu0; /* diskheight is height of disk surface */
  double density_at_surface;
  int ngrid;
  double *zht, *temp, *den, *tau;
};

/* structure containing opacity information */
struct disk_opacity {
  /* opacity_d -> Rosseland,  opacity_V -> Planck */
  //  double opacity_d, opacity_V;
  double Rosseland, Planck;
  double ratio, ratio_p;  /* ratio = Planck/Rosseland, ratio_p = */
  double abs_frac, sca_frac, g_param;
};

/* structure containing calculated temperature data */
struct local_perturbation {
  int xgrid, ygrid;
  double xmax, phimax;
  double *xpos, *phipos;
  double **flux_disk, **surface; //just 2D arrays now
};

/* structure containing variables that can be calculated before flux loop */
struct var_struct {
  int *non_diffuse;
  double *normx, *normy, *normz;
  double *midx, *midy, *midz;
  double *F_irr;
  double *omeg;
  double *d0_starx, *d0_stary, *d0_starz; 
  double *c1p, *c2p, *c3p;
  double *mu;
  /*double *mu_std, *omeg_box;
  double *normx_box, *normy_box, *normz_box;
  double *midx_box, *midy_box, *midz_box;*/
};


		      
void surface_range(struct disk_parameters *DP, 
                 struct disk_opacity *opacity, struct var_struct *VS,
                 double *full_surf, double *full_dens, double *rpos);
void fix_ghost_cells();
void freeVS(struct var_struct *VS);                 		      

double cartesian2cylindrical(double rx, double ry, double *r,double *phi);
double cylindrical2cartesian(double r, double phi, double *rx,double *ry);
double array_interp(double x, double y, double r0, 
                    double dr, double dphi, double *arr);
double array_interp2(double x, double y, double r0, double *rpos, 
                    double dphi, double *arr, double *arr2, double *val2);  
double array_interp4(double x, double y, double r0, double *rpos, 
                    double dphi, double *arr, double *arr2, double *arr3, float *arr4,
                    double *val2, double *val3, float *val4);  
double B_flux_azi(struct var_struct *VS, double *full_H, double *full_dens,
                     int i, int j, int ii, int jj, int i_shift, int ll,
                     double *nu_omega, int ygrid, double ox, double oy, double oz,
                     double opac, double g_param, double r, double dr, double r0,
                     double *tau_array, int *flux_method);


void calc_flux_azi5(struct disk_parameters *DP,struct disk_opacity *opacity,
                   struct var_struct *VS, double *full_dens, double *full_H,
                   double *full_surf, float *full_cs, double *flux_azi, double *dens_azi,
                   int opac_i, double *rpos);	
double B_flux_azi5(struct var_struct *VS, double *full_H, double *full_dens,
                     int i, int j, int ii, int jj, int i_shift, int ll,
                     double *nu_omega, int ygrid, double ox, double oy, double oz,
                     struct disk_opacity *opacity, double r, double *rpos, double r0,
                     double *tau_return, int *flux_method, double tau_min,
                     double ph_nx, double ph_ny, double ph_nz, double zph);      
void viscous_accretion_flux(double *flux_sub, struct disk_parameters *DP, 
                            struct disk_opacity *opacity);

// used for extrapolating surface 
#define NFIT 8
double fit_alpha(int npts, double rvals[], double hvals[]);

// useful functions 
double F_incident(double T_star, double r_star, double a_AU);
double Hill_radius(double M_star, double M_planet, double a_au);
double Incident_Slope(double mu0);

// disk geometry 
//void calcSurface(double H0, double rhill, double alpha, double distance);
double incident_angle(double rmid, double zmid,
		      double dzdr, double dzdphi, double R_star);
double min_incident_angle(double dtot, double R_star);

double get_photosphere(double sigma, double H, double opac, double depth, float limit);

double get_tau_min(int ii, int jj,double *full_dens,double *full_H,double opac,
                   double *rpos, double r0, struct var_struct *VS);

void calc_surf_twotemp(struct disk_parameters *DP, struct disk_opacity *opacity,
                   double *dens_azi, double *H_azi, double *surf_azi,int opac_i,
                   int size_x,int size_y,double avg_rho_in, double avg_alpha_rho, double *rpos);
void calc_surf_lookup(struct disk_parameters *DP, struct disk_opacity *opacity,
                   double *dens_azi, double *H_azi, double *surf_azi,int opac_i,
                   int size_x,int size_y,double avg_rho_in, double *rpos);
void output_density_struct(struct disk_parameters *DP, struct disk_opacity *opacity,
                    double *dens_azi, double *H_azi, double *surf_azi,int opac_i,
                    int size_x,int size_y, double *rpos);
void surface_smoothing(double *surf_azi);

double calc_upper_incomp_gamma(double s, double x);

double integrate_gauss(double x1, double x2, int printint);

double moment_solidangle_strip(double r0, double z0,
			       double rin, double zin,
			       double rout, double zout);

void calc_flux_azi(struct disk_parameters *DP,struct disk_opacity *opacity,
                   struct var_struct *VS, double *full_dens, double *full_H,
                   double *full_surf,double *flux_azi, double *dens_azi,
                   int opac_i);	       


double B_inner_disk_HSL(int j, int j_shift, double *full_surf, double *rpos,
                       double *nu_omega_tot, struct disk_parameters *DP, 
                       struct disk_opacity *opacity, struct var_struct *VS, int cells_in);
double B_outer_disk_HSL(int j, int j_shift, double *full_surf, double *rpos,
                      double *nu_omega_tot, struct disk_parameters *DP, 
                      struct disk_opacity *opacity, struct var_struct *VS);


double B_unpert_los(double tau_los, double mu, double F_irr, 
		    struct disk_opacity *opacity);
double B_unpert_diffuse(double mu, double F_irr, 
		    struct disk_opacity *opacity);
		    
  
typedef struct _vector{
  double x,y,z;
} vector;

typedef struct _point{
  double x,y,z;
} point;

#define SQUARE 4
double dotproduct(vector v1, vector v2);
#define cosangle(v1,v2) ( dotproduct(v1,v2)/sqrt(dotproduct(v1,v1)*dotproduct(v2,v2)) )
vector unitvector(vector v);
vector crossproduct(vector v1, vector v2);
vector vectorbtwnpts(point p1, point p2);
point cylin2cart(double x, double phi, double z, double a_au);
point pointat(double x, double y, double z);		    

int gsl_fit_wlinear (const double *x, const size_t xstride,
                 const double *w, const size_t wstride,
                 const double *y, const size_t ystride,
                 const size_t n,
                 double *c0, double *c1,
                 double *cov_00, double *cov_01, double *cov_11,
                 double *chisq);

// for printf-ing comments 
//#define MPI_printf(RANK,MY_STRING) ( if(RANK==0) printf(MY_STRING); )

ex void init_stockholm_HSL(double *full_dens, double *rpos);
//ex void StockholmBoundary_HSL_cpu(real);

