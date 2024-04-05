/*accuracy: 8th-order space, 2nd-order time*/
/*avoid low-order FD at near-boundary grids, which can introduce dispersion*/

#include "funclist.h"

void step_forward(int nxpad, int nypad, float dx, float dy, float dt, int fdo, float *hc,
                  float **inv_rho_half_x, float **inv_rho_half_y, float **rhopad, float **vppad,
                  float *K_x, float *a_x, float *b_x, float *K_y, float *a_y, float *b_y,
                  float *K_x_half, float *a_x_half, float *b_x_half, float *K_y_half, float *a_y_half, float *b_y_half, 
                  float **psi_dp_dx, float **psi_dp_dy, float **psi_dvx_dx, float **psi_dvy_dy, float **p, float **v_x, float **v_y)

{
  int i,j;
  float value_dp_dx, value_dp_dy, value_dvx_dx, value_dvy_dy;
  
  

  //********update particle velocities***************
  for ( j=fdo; j<nypad+fdo; j++ ) {
	for ( i=fdo; i<nxpad+fdo; i++ ) {

		value_dp_dx = (hc[1]* ( p[j][i+1]-p[j][i] )+
		               hc[2]* ( p[j][i+2]-p[j][i-1] )+
			       hc[3]* ( p[j][i+3]-p[j][i-2] )+
			       hc[4]* ( p[j][i+4]-p[j][i-3] ))/dx;

		value_dp_dy = (hc[1]* ( p[j+1][i]-p[j][i] )+ 
			       hc[2]* ( p[j+2][i]-p[j-1][i] )+
			       hc[3]* ( p[j+3][i]-p[j-2][i] )+
			       hc[4]* ( p[j+4][i]-p[j-3][i] ))/dy;
		
		//left-right CPML
		psi_dp_dx[j-fdo][i-fdo] = b_x_half[i-fdo] * psi_dp_dx[j-fdo][i-fdo] + a_x_half[i-fdo] * value_dp_dx;
		value_dp_dx = value_dp_dx / K_x_half[i-fdo] + psi_dp_dx[j-fdo][i-fdo];
	
		//top-bottom CPML	
		psi_dp_dy[j-fdo][i-fdo] = b_y_half[j-fdo] * psi_dp_dy[j-fdo][i-fdo] + a_y_half[j-fdo] * value_dp_dy;
		value_dp_dy = value_dp_dy / K_y_half[j-fdo] + psi_dp_dy[j-fdo][i-fdo];
		
		v_x[j][i] -= dt*value_dp_dx *inv_rho_half_x[j-fdo][i-fdo];
		v_y[j][i] -= dt*value_dp_dy *inv_rho_half_y[j-fdo][i-fdo];

	}
  }

  //******************pressure update****************************
  for ( j=fdo; j<nypad+fdo; j++ ) {
	for ( i=fdo; i<nxpad+fdo; i++ ) {
		value_dvx_dx = ( hc[1]* ( v_x[j][i]  -v_x[j][i-1] )+
		                 hc[2]* ( v_x[j][i+1]-v_x[j][i-2] )+
		                 hc[3]* ( v_x[j][i+2]-v_x[j][i-3] )+
		                 hc[4]* ( v_x[j][i+3]-v_x[j][i-4] )) /dx;
		
		value_dvy_dy = ( hc[1]* ( v_y[j][i]  -v_y[j-1][i] )+
			         hc[2]* ( v_y[j+1][i]-v_y[j-2][i] )+
			         hc[3]* ( v_y[j+2][i]-v_y[j-3][i] )+
			         hc[4]* ( v_y[j+3][i]-v_y[j-4][i] )) /dy;

		//left-right CPML
		psi_dvx_dx[j-fdo][i-fdo] = b_x[i-fdo] * psi_dvx_dx[j-fdo][i-fdo] + a_x[i-fdo] * value_dvx_dx;
		value_dvx_dx = value_dvx_dx / K_x[i-fdo] + psi_dvx_dx[j-fdo][i-fdo];

		//top-bottom CPML
		psi_dvy_dy[j-fdo][i-fdo] = b_y[j-fdo] * psi_dvy_dy[j-fdo][i-fdo] + a_y[j-fdo] * value_dvy_dy;
		value_dvy_dy = value_dvy_dy / K_y[j-fdo] + psi_dvy_dy[j-fdo][i-fdo];
		
		p[j][i] -= dt* (value_dvx_dx + value_dvy_dy) * rhopad[j-fdo][i-fdo]*vppad[j-fdo][i-fdo]*vppad[j-fdo][i-fdo];
	}
  }
  


}






