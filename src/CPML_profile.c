
#include "funclist.h"


void CPML_profile(float * d_x, float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
            float * d_x_half, float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
            float * d_z, float * K_z, float * alpha_prime_z, float * a_z, float * b_z, 
            float * d_z_half, float * K_z_half, float * alpha_prime_z_half, float * a_z_half, float * b_z_half, 
            int nxpad, int nzpad, int nb, float dx, float dz, float VPPML, float dt, float f0)
{

      /* local variables */
      int i;
      
      float NPOWER=2.0, K_MAX_CPML=1.0;
      const float alpha_max_PML = 2.0 * PI * (f0/2.0); /* from festa and Vilotte */

      float thickness_PML_x, thickness_PML_z, xoriginleft, xoriginright, zoriginbottom, zorigintop;
      float Rcoef, d0_x, d0_z, xval, zval, abscissa_in_PML, abscissa_normalized;


      /* thickness of the PML layer in meters */
      thickness_PML_x = (float)nb*dx;
      thickness_PML_z = (float)nb*dz;

      /* reflection coefficient (INRIA report section 6.1) */
      Rcoef = 0.001;

      /* compute d0 from INRIA report section 6.1 */
      d0_x = - (NPOWER + 1) * VPPML * log(Rcoef) / (2.0 * thickness_PML_x);
      d0_z = - (NPOWER + 1) * VPPML * log(Rcoef) / (2.0 * thickness_PML_z);
	
      /* damping in the X direction */

      xoriginleft = thickness_PML_x;
      xoriginright = (nxpad-1) * dx - thickness_PML_x;

      for (i=0;i<nxpad;i++){       
          
            xval = dx * (float)(i);   //xval = dx * (float)(i-1);
            
            /* left boundary */
            abscissa_in_PML = xoriginleft - xval;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x[i] = d0_x * pow(abscissa_normalized,NPOWER);

            /* this taken from Gedney page 8.2 */
               K_x[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_x[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = xoriginleft - (xval + dx/2.0);

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x_half[i] = d0_x * pow(abscissa_normalized,NPOWER);

               /* this taken from Gedney page 8.2 */
               K_x_half[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_x_half[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

            /* right boundary */                                
            abscissa_in_PML = xval - xoriginright;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x[i] = d0_x * pow(abscissa_normalized,NPOWER);

            /* this taken from Gedney page 8.2 */
               K_x[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_x[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = xval + dx/2.0 - xoriginright;

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_x;
               d_x_half[i] = d0_x * pow(abscissa_normalized,NPOWER);

               /* this taken from Gedney page 8.2 */
               K_x_half[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_x_half[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

            /* just in case, for -5 at the end */
            if(alpha_prime_x[i] < 0.0){ alpha_prime_x[i] = 0.0;}
            if(alpha_prime_x_half[i] < 0.0) {alpha_prime_x_half[i] = 0.0;}

            b_x[i] = exp(- (d_x[i] / K_x[i] + alpha_prime_x[i]) * dt);
            b_x_half[i] = exp(- (d_x_half[i] / K_x_half[i] + alpha_prime_x_half[i]) * dt);

            /* avoid division by zero outside the PML */
            if(fabs(d_x[i]) > 1.0e-6){ a_x[i] = d_x[i] * (b_x[i] - 1.0) / (K_x[i] * (d_x[i] + K_x[i] * alpha_prime_x[i]));}
            if(fabs(d_x_half[i]) > 1.0e-6){ a_x_half[i] = d_x_half[i] * (b_x_half[i] - 1.0) / (K_x_half[i] * (d_x_half[i] + K_x_half[i] * alpha_prime_x_half[i]));}

      }   



	
      /* damping in the Z direction */

      zorigintop = thickness_PML_z;
      zoriginbottom = (nzpad-1)*dz - thickness_PML_z;

      for (i=0;i<nzpad;i++){
            zval = dz * (float)(i);

          /* top boundary */
            abscissa_in_PML = zorigintop - zval;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_z;
               d_z[i] = d0_z * pow(abscissa_normalized,NPOWER);

            /* this taken from Gedney page 8.2 */
               K_z[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_z[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = zorigintop - (zval + dz/2.0);

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_z;
               d_z_half[i] = d0_z * pow(abscissa_normalized,NPOWER);

               /* this taken from Gedney page 8.2 */
               K_z_half[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_z_half[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

      /* bottom boundary */      
            abscissa_in_PML = zval - zoriginbottom;
      
            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_z;
               d_z[i] = d0_z * pow(abscissa_normalized,NPOWER);

            /* this taken from Gedney page 8.2 */
               K_z[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_z[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

            /* define damping profile at half the grid points */
            abscissa_in_PML = zval + dz/2.0 - zoriginbottom;

            if(abscissa_in_PML >= 0.0){
               abscissa_normalized = abscissa_in_PML / thickness_PML_z;
               d_z_half[i] = d0_z * pow(abscissa_normalized,NPOWER);

               /* this taken from Gedney page 8.2 */
               K_z_half[i] = 1.0 + (K_MAX_CPML - 1.0) * pow(abscissa_normalized,NPOWER);
               alpha_prime_z_half[i] = alpha_max_PML * (1.0 - abscissa_normalized); // + 0.1 * alpha_max_PML
            }

          b_z[i] = exp(- (d_z[i] / K_z[i] + alpha_prime_z[i]) * dt);
          b_z_half[i] = exp(- (d_z_half[i] / K_z_half[i] + alpha_prime_z_half[i]) * dt);

          /* avoid division by zero outside the PML */
          if(fabs(d_z[i]) > 1.0e-6){ a_z[i] = d_z[i] * (b_z[i] - 1.0) / (K_z[i] * (d_z[i] + K_z[i] * alpha_prime_z[i]));}
          if(fabs(d_z_half[i]) > 1.0e-6){ a_z_half[i] = d_z_half[i] * (b_z_half[i] - 1.0) / (K_z_half[i] * (d_z_half[i] + K_z_half[i] * alpha_prime_z_half[i]));}
      
      }  

}
