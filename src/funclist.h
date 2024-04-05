#include <rsf.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>


#define PI 3.141592653589793238462643



void step_forward(int nxpad, int nypad, float dx, float dy, float dt, int fdo, float *hc,
                  float **inv_rho_half_x, float **inv_rho_half_y, float **rhopad, float **vppad,
                  float *K_x, float *a_x, float *b_x, float *K_y, float *a_y, float *b_y,
                  float *K_x_half, float *a_x_half, float *b_x_half, float *K_y_half, float *a_y_half, float *b_y_half, 
                  float **psi_dp_dx, float **psi_dp_dy, float **psi_dvx_dx, float **psi_dvy_dy, float **p, float **v_x, float **v_y);
                  
void bndr_rw(bool read, float **p, char *tmp_srcwfd_dir, int it, int isrc, int fdo, int nb, int nx, int nz);

void outsnap( int snaptype, float **p, float **vx, float **vz, int isrc, int it, int nxpad, int nzpad, int nb, int fdo, char *snapdir, char ext[50] );

void outsnapimg(float **num, float **den, int it, int MYID, int nx, int nz, char *snapimg_dir);

void cross_correlation( float **num, float **den, float **sp, float **gp, int nb, int fdo, int nx, int nz );

void matrix_transp_2D(float **mat_ex, float **mat_in, int nx, int nz, int io);

void pad_bound ( int nx, int nz, int nxpad, int nzpad, int nb, float **raw_mod, float **pad_mod );

void CPML_profile(float * d_x, float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
            float * d_x_half, float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
            float * d_z, float * K_z, float * alpha_prime_z, float * a_z, float * b_z, 
            float * d_z_half, float * K_z_half, float * alpha_prime_z_half, float * a_z_half, float * b_z_half, 
            int nxpad, int nzpad, int nb, float dx, float dz, float VPPML, float dt, float f0);

float coef(float *hc, int maxerror, int fdorder, int MYID);

int **source(int *nsrc, char *SOURCE_FILE, int MYID, float dx, float dz, int nx, int nz, int nb);

int **receiver(int *ntr, char *REC_FILE, int MYID, float dx, float dz, int nx, int nz, int nb);







