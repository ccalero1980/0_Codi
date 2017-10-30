#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nr.h"
#include "nrutil.h"
#include "nrutil.c"
//#include "vectors.h"
//#include "vectors.c"
#include "tridag.c"
#include "cyclic.c"


/* ######################################################################### */

void  Torque_Bfield( vec_s  r1,
           vec_s  r2,
           vec_s  m2,
           vec_s  *f_Btorque,
           param_s params,
           vec_s Bfield)
/* ######################################################################### */
{
  float r21x, r21z, Torque; 
  vec_s res;  

  Torque   =  (Bfield.x*m2.z - Bfield.z*m2.x); // Suposo que el camp magnetic es al pla XZ  

  r21x =  (r2.x - r1.x);
  r21z =  (r2.z - r1.z);

  f_Btorque[0].x = r21z * Torque/(params.b_l * params.b_l);
  f_Btorque[0].z = -r21x * Torque/(params.b_l * params.b_l);

  f_Btorque[0].x /= params.mass;
  f_Btorque[0].z /= params.mass;
  
}

/* ######################################################################### */

void  Torque_Bfield_2filaments( vec_s  r1,
           vec_s  r2,
           vec_s  m2,
           vec_s  *f_Btorque,
           param_s params,
           vec_s Bfield)
/* ######################################################################### */
{
  float r21x, r21z, Torque; 
  vec_s res;  

  Torque   =  (Bfield.x*m2.z - Bfield.z*m2.x); // Suposo que el camp magnetic es al pla XZ  

  r21x = - (r2.x - r1.x)/2.0;
  r21z =  (r2.z - r1.z)/2.0;

  f_Btorque[0].x = 0.5*r21z*Torque/((params.b_l/2.)*(params.b_l/2.));
  f_Btorque[0].z = 0.5*r21x*Torque/((params.b_l/2.)*(params.b_l/2.));

  f_Btorque[0].x /= params.mass1;
  f_Btorque[0].z /= params.mass1;
  
}

/* ######################################################################### */

vec_s  Binduction( vec_s  r0f,
           vec_s  m2, 
           param_s params)
/* ######################################################################### */
{
  vec_s v1, v2, Bin;
  float mu0, sc1, d2, d, d3, d5;

  mu0 = params.mu0;

  d2 = r0f.x*r0f.x + r0f.y*r0f.y + r0f.z*r0f.z;
  d = sqrt(d2);
  d3 = d2*d;
  d5 = d2*d2*d;

  sc1 = dotProduct(r0f,m2);
  v1 = scalarProduct(r0f, 3.*sc1/d5);
  v2 = scalarProduct(m2, -1./d3);
  Bin = vecAddition(v1, v2);
  Bin = scalarProduct(Bin, mu0/(4.*PI));

  return(Bin);

}

/* ######################################################################### */

void  dipolar_forces( vec_s  r0f,
           vec_s  m1,
           vec_s  m2,
           vec_s  *f_dip,
           param_s params)
/* ######################################################################### */
{

    int   n_mon, i;
    float mu0, d2, d, d5, sc1;
    vec_s res, v1, v2;
    res.x = 0.0;
    res.y = 0.0;
    res.z = 0.0;
    
    mu0 = params.mu0;

    
    d2 = r0f.x*r0f.x + r0f.y*r0f.y + r0f.z*r0f.z;
    d5 = sqrt(d2)*d2*d2;
    
    v1 = crossProduct(r0f, m1);
    v2 = crossProduct(v1, m2);    
    res = vecAddition(res, v2);
    
    v1 = crossProduct(r0f, m2);
    v2 = crossProduct(v1, m1);
    res = vecAddition(res, v2);

    sc1 = dotProduct(m1,m2);
    sc1 = -2*sc1;
    v2 = scalarProduct(r0f, sc1);
    res = vecAddition(res, v2);
    
    v1 = crossProduct(r0f, m1);
    v2 = crossProduct(r0f, m2);   
    sc1 = dotProduct(v1,v2);
    sc1 = 5*sc1/d2;
    v2 = scalarProduct(r0f, sc1);
    res = vecAddition(res, v2);
    
    f_dip[0].x = res.x*3.*mu0/(4*PI*d5);
    f_dip[0].y = res.y*3.*mu0/(4*PI*d5);
    f_dip[0].z = res.z*3.*mu0/(4*PI*d5);
  
}

/* ######################################################################### */

void  dipolar_forces_2filaments( vec_s  r0f,
           vec_s  m2,
           vec_s  *f_dip,
           param_s params,
           vec_s Bfield)
/* ######################################################################### */
{

    int   n_mon, i;
    float mu0, d2, d, d5, sc1;
    vec_s m1, res, v1, v2;
    res.x = 0.0;
    res.y = 0.0;
    res.z = 0.0;
    
    mu0 = params.mu0;

    m1.x = params.chi*Bfield.x;
    m1.y = params.chi*Bfield.y;
    m1.z = params.chi*Bfield.z;
    
    d2 = r0f.x*r0f.x + r0f.y*r0f.y + r0f.z*r0f.z;
    d5 = sqrt(d2)*d2*d2;
    
    v1 = crossProduct(r0f, m1);
    v2 = crossProduct(v1, m2);    
    res = vecAddition(res, v2);
    
    v1 = crossProduct(r0f, m2);
    v2 = crossProduct(v1, m1);
    res = vecAddition(res, v2);

    sc1 = dotProduct(m1,m2);
    sc1 = -2*sc1;
    v2 = scalarProduct(r0f, sc1);
    res = vecAddition(res, v2);
    
    v1 = crossProduct(r0f, m1);
    v2 = crossProduct(r0f, m2);   
    sc1 = dotProduct(v1,v2);
    sc1 = 5*sc1/d2;
    v2 = scalarProduct(r0f, sc1);
    res = vecAddition(res, v2);
    
    f_dip[0].x += res.x*3.*mu0/(4*PI*d5);
    f_dip[0].y += res.y*3.*mu0/(4*PI*d5);
    f_dip[0].z += res.z*3.*mu0/(4*PI*d5);
  
}

/* ######################################################################### */

void  dipolar_forces_2filaments_induction( vec_s  r0f,
           vec_s  m1,
           vec_s  m2,
           vec_s  *f_dip,
           param_s params,
           vec_s Bfield)
/* ######################################################################### */
{

    int   n_mon, i;
    float mu0, d2, d5, sc1;
    vec_s res, v1, v2;
    res.x = 0.0;
    res.y = 0.0;
    res.z = 0.0;
    mu0 = params.mu0;

    d2 = r0f.x*r0f.x + r0f.y*r0f.y + r0f.z*r0f.z;
    d5 = sqrt(d2)*d2*d2;
    
    v1 = crossProduct(r0f, m1);
    v2 = crossProduct(v1, m2);    
    res = vecAddition(res, v2);
    
    v1 = crossProduct(r0f, m2);
    v2 = crossProduct(v1, m1);
    res = vecAddition(res, v2);

    sc1 = dotProduct(m1,m2);
    sc1 = -2*sc1;
    v2 = scalarProduct(r0f, sc1);
    res = vecAddition(res, v2);
    
    v1 = crossProduct(r0f, m1);
    v2 = crossProduct(r0f, m2);   
    sc1 = dotProduct(v1,v2);
    sc1 = 5*sc1/d2;
    v2 = scalarProduct(r0f, sc1);
    res = vecAddition(res, v2);
    
    f_dip[0].x += res.x*3.*mu0/(4*PI*d5);
    f_dip[0].y += res.y*3.*mu0/(4*PI*d5);
    f_dip[0].z += res.z*3.*mu0/(4*PI*d5);
  
}




/* ######################################################################### */

void  steric_force( vec_s  r0f, 
                    float  a,
                    vec_s *f_steric)
/* ######################################################################### */
{
    float d2, inv2, inv4, inv8, inv14;
    vec_s res;
    d2 = dotProduct(r0f,r0f);
    inv2 = a*a/d2;
    inv4 = inv2*inv2;
    inv8 = inv4*inv4;
    inv14 = inv8*inv4*inv2;
    res = scalarProduct(r0f, inv14);

    if(inv2>1.0){
      f_steric[0].x += res.x; 
      f_steric[0].y += res.y;
      f_steric[0].z += res.z;
    }
    else{
      f_steric[0].x += 0.0; 
      f_steric[0].y += 0.0;
      f_steric[0].z += 0.0;      
    }
    //printf("dins steric %le\t%le\t%le\t%le\t%le\t%le\n", d2, a, inv8, f_steric[0].x, f_steric[0].y,f_steric[0].z);

 
  }

/* ######################################################################### */

void  steric_force2( vec_s  r0f, 
                    float  a,
                    vec_s *f_steric)
/* ######################################################################### */
{
    float d2, inv2, sc1, inv6, inv12;
    a = a/pow(2.,1./6.);
    vec_s res;
    d2 = dotProduct(r0f,r0f);
    inv2 = a*a/d2;
    inv6 = inv2*inv2*inv2;
    inv12 = inv6*inv6;
    sc1 = (12.*inv12 - 6.*inv6)/d2;
    res = scalarProduct(r0f, sc1);

    if(inv6>0.5){
      f_steric[0].x += res.x; 
      f_steric[0].y += res.y;
      f_steric[0].z += res.z;
    }
    else{
      f_steric[0].x += 0.0; 
      f_steric[0].y += 0.0;
      f_steric[0].z += 0.0;      
    }
    //printf("dins steric %le\t%le\t%le\t%le\t%le\t%le\n", d2, a, inv8, f_steric[0].x, f_steric[0].y,f_steric[0].z);

 
  }
/* ######################################################################### */

float  bending_forces( vec_s  *rij1,
		       float  *mags1,
		       vec_s  *f,
		       param_s params    )
/* ######################################################################### */
{


    int    n_mon, i;
   
    float  b_l, k, phi_0, u, m_inv;
 
    vec_s  f_di[3];
  

    n_mon = params.n_mon;
    b_l   = params.b_l;
    k     = params.k;
    phi_0 = params.phi_0;
    u     = 0.0;
    m_inv = 1.0/params.mass;
    for( i = 0; i < n_mon - 2; i++)
    {
       u += angle_force( rij1[i + 1], 
                         rij1[i+2], 
                         mags1[i+1], 
                         mags1[i+2],
                         k,
                         phi_0,
                         f_di);

       f[i].x += f_di[0].x*m_inv;
       f[i].y += f_di[0].y*m_inv;
       f[i].z += f_di[0].z*m_inv;
       f[i+1].x += f_di[1].x*m_inv;
       f[i+1].y += f_di[1].y*m_inv;
       f[i+1].z += f_di[1].z*m_inv;
       f[i+2].x += f_di[2].x*m_inv;
       f[i+2].y += f_di[2].y*m_inv;
       f[i+2].z += f_di[2].z*m_inv;
    }
    return( u );  
}


/* ######################################################################### */

float  bending_forces_2filaments( vec_s  *rij1,
           float  *mags1,
           vec_s  *f,
           param_s params    )
/* ######################################################################### */
{


    int    n_mon1, i;
   
    float  b_l, k, phi_0, u, m_inv;
 
    vec_s  f_di[3];
  

    n_mon1 = params.n_mon1;
    b_l   = params.b_l;
    k     = params.k1;
    phi_0 = params.phi_0;
    u     = 0.0;
    m_inv = 1.0/params.mass1;
    for( i = 0; i < n_mon1 - 2; i++)
    {
       u += angle_force( rij1[i + 1], 
                         rij1[i+2], 
                         mags1[i+1], 
                         mags1[i+2],
                         k,
                         phi_0,
                         f_di);

       f[i].x += f_di[0].x*m_inv;
       f[i].y += f_di[0].y*m_inv;
       f[i].z += f_di[0].z*m_inv;
       f[i+1].x += f_di[1].x*m_inv;
       f[i+1].y += f_di[1].y*m_inv;
       f[i+1].z += f_di[1].z*m_inv;
       f[i+2].x += f_di[2].x*m_inv;
       f[i+2].y += f_di[2].y*m_inv;
       f[i+2].z += f_di[2].z*m_inv;
    }
    return( u );  
}


/* ######################################################################### */

float  bending_forces_circ( vec_s  *rij1,
           float  *mags1,
           vec_s  *f,
           param_s params    )
/* ######################################################################### */
{


    int    n_mon, i;
   
    float  b_l, k, phi_0, u, m_inv;
 
    vec_s  f_di[3];
  

    n_mon = params.n_mon;
    b_l   = params.b_l;
    k     = params.k;
    phi_0 = params.phi_0;
    u     = 0.0;
    m_inv = 1.0/params.mass;
    for( i = 0; i < n_mon - 2; i++)
    {
       u += angle_force( rij1[i + 1], 
                         rij1[i+2], 
                         mags1[i+1], 
                         mags1[i+2],
                         k,
                         phi_0,
                         f_di);

       f[i].x += f_di[0].x*m_inv;
       f[i].y += f_di[0].y*m_inv;
       f[i].z += f_di[0].z*m_inv;
       f[i+1].x += f_di[1].x*m_inv;
       f[i+1].y += f_di[1].y*m_inv;
       f[i+1].z += f_di[1].z*m_inv;
       f[i+2].x += f_di[2].x*m_inv;
       f[i+2].y += f_di[2].y*m_inv;
       f[i+2].z += f_di[2].z*m_inv;
    }
       u += angle_force( rij1[n_mon - 1], 
                         rij1[n_mon], 
                         mags1[n_mon-1], 
                         mags1[n_mon],
                         k,
                         phi_0,
                         f_di);


       f[n_mon-2].x += f_di[0].x*m_inv;
       f[n_mon-2].y += f_di[0].y*m_inv;
       f[n_mon-2].z += f_di[0].z*m_inv;
       f[n_mon-1].x += f_di[1].x*m_inv;
       f[n_mon-1].y += f_di[1].y*m_inv;
       f[n_mon-1].z += f_di[1].z*m_inv;
       f[0].x += f_di[2].x*m_inv;
       f[0].y += f_di[2].y*m_inv;
       f[0].z += f_di[2].z*m_inv;

       u += angle_force( rij1[n_mon], 
                         rij1[1], 
                         mags1[n_mon], 
                         mags1[1],
                         k,
                         phi_0,
                         f_di);


       f[n_mon-1].x += f_di[0].x*m_inv;
       f[n_mon-1].y += f_di[0].y*m_inv;
       f[n_mon-1].z += f_di[0].z*m_inv;
       f[0].x += f_di[1].x*m_inv;
       f[0].y += f_di[1].y*m_inv;
       f[0].z += f_di[1].z*m_inv;
       f[1].x += f_di[2].x*m_inv;
       f[1].y += f_di[2].y*m_inv;
       f[1].z += f_di[2].z*m_inv;

    return( u );  
}

/* ########################################################################## */
float   angle_force( vec_s   r1,
                    vec_s   r2, 
                    float   r12,
                    float   r23,
                    float   k,
                    float   phi_0,
                    vec_s  *f )
/* ########################################################################## */
{
   int   i;

   float dprod, fac, u;
   float r12_3, r23_3, cphi, sphi, arg;

   dprod = r1.x*r2.x + r1.y*r2.y + r1.z*r2.z;
   cphi   = dprod/(r12*r23);

   u      = k*(1.0 - cphi);
   fac    = k;

   /*     sphi   = (1.0 - sqrt(cphi))/2.0; */
   /*     u      = k*(sphi)/2.0; */
   /*     fac    = k/4.0; */
 
   r12_3  = r12*r12*r12;
   r23_3  = r23*r23*r23;


   f[0].x =  r2.x/(r12*r23) - r1.x*dprod/(r12_3*r23);
   f[0].y =  r2.y/(r12*r23) - r1.y*dprod/(r12_3*r23);
   f[0].z =  r2.z/(r12*r23) - r1.z*dprod/(r12_3*r23);
   f[1].x = -r2.x*dprod/(r12*r23_3) + 
            (r1.x - r2.x)/(r12*r12) +
             r1.x*dprod/(r12_3*r23);
   f[1].y = -r2.y*dprod/(r12*r23_3) + 
            (r1.y - r2.y)/(r12*r12) +
             r1.y*dprod/(r12_3*r23);
   f[1].z = -r2.z*dprod/(r12*r23_3) + 
            (r1.z - r2.z)/(r12*r12) +
             r1.z*dprod/(r12_3*r23);
   f[2].x  = r2.x*dprod/(r12*r23_3) - r1.x/(r12*r23); 
   f[2].y  = r2.y*dprod/(r12*r23_3) - r1.y/(r12*r23); 
   f[2].z  = r2.z*dprod/(r12*r23_3) - r1.z/(r12*r23);
   for( i = 0; i < 3; i++ )
   {
      f[i].x *= fac;
      f[i].y *= fac;
      f[i].z *= fac;
   }
   return( u );
  
}
/* ########################################################################## */
void     constrain_positions( vec_s  *r1,
                              vec_s  *rij1,
                              vec_s  *v1,
                              vec_s  *fc_2,
                              param_s params )
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans;

    int   n_mon, n_c, i;
    float rijx, rijy, rijz, mag;
    float rhs_old, b_l, b_l_sq, max_error, error, tol;
    float dx, dy, dz, dt, disp1;

    n_mon   = params.n_mon;
    n_c     = n_mon - 1;
    tol     = params.tolerance;
    dt      = params.dt;
    b_l     = params.b_l;
    b_l_sq  = b_l*b_l;
    disp1    = 1.0/(1.0 + 0.5*params.gamma*params.dt);

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_2[i].x = 0.0;
        fc_2[i].y = 0.0;
        fc_2[i].z = 0.0;
    }
    ans[0]     = 0.0;
    ans[n_mon] = 0.0;
    rij1[0].x   = 0.0;
    rij1[0].y   = 0.0;
    rij1[0].z   = 0.0;
    rij1[n_c + 1].x = 0.0;
    rij1[n_c + 1].y = 0.0;
    rij1[n_c + 1].z = 0.0;

    rijx  = r1[0].x - r1[1].x;
    rijy  = r1[0].y - r1[1].y;
    rijz  = r1[0].z - r1[1].z;
    mag   = rijx*rijx + rijy*rijy + rijz*rijz; 

    a[1]   = 0.0;
    b[1]   = rij1[1].x*rijx;
    b[1]  += rij1[1].y*rijy;
    b[1]  += rij1[1].z*rijz;
    b[1]  *= 4.0;
    c[1]   = rij1[2].x*rijx;
    c[1]  += rij1[2].y*rijy;
    c[1]  += rij1[2].z*rijz;
    c[1]  *= -2.0;
    rhs[1] = b_l_sq - mag;

    rijx    = r1[n_c - 1].x - r1[n_c].x;
    rijy    = r1[n_c - 1].y - r1[n_c].y;
    rijz    = r1[n_c - 1].z - r1[n_c].z;
    mag     = rijx*rijx + rijy*rijy + rijz*rijz; 
    a[n_c]  = rij1[n_c - 1].x*rijx;
    a[n_c] += rij1[n_c - 1].y*rijy;
    a[n_c] += rij1[n_c - 1].z*rijz;
    a[n_c] *= -2.0;
    b[n_c]  = rij1[n_c].x*rijx;
    b[n_c] += rij1[n_c].y*rijy;
    b[n_c] += rij1[n_c].z*rijz;
    b[n_c] *= 4.0;
    c[n_c]  = 0.0;
    rhs[n_c]= b_l_sq - mag;
    for( i = 2; i < n_c; i++ )
    {
        rijx   = r1[i-1].x - r1[i].x;
        rijy   = r1[i-1].y - r1[i].y;
        rijz   = r1[i-1].z - r1[i].z;
        mag    = rijx*rijx + rijy*rijy + rijz*rijz; 
        a[i]   = -2.0*(rijx*rij1[i-1].x + rijy*rij1[i-1].y + rijz*rij1[i-1].z);
        b[i]   =  4.0*(rijx*rij1[i].x + rijy*rij1[i].y + rijz*rij1[i].z);
        c[i]   = -2.0*(rijx*rij1[i+1].x + rijy*rij1[i+1].y + rijz*rij1[i+1].z);
        rhs[i] = b_l_sq - mag;
    }
    do
    {  
       tridag( a, b, c, rhs, ans, n_c );
       max_error = 0.0;
       for( i = 1; i <= n_c; i++ )
       {
          rijx    = r1[i-1].x + rij1[i].x*ans[i] - rij1[i-1].x*ans[i-1] -
                   (r1[i].x + rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i]);
          rijy    = r1[i-1].y + rij1[i].y*ans[i] - rij1[i-1].y*ans[i-1] -
                   (r1[i].y + rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i]);
          rijz    = r1[i-1].z + rij1[i].z*ans[i] - rij1[i-1].z*ans[i-1] -
                   (r1[i].z + rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i]);
          mag     = rijx*rijx + rijy*rijy + rijz*rijz;
          error   = (sqrt(mag) - b_l)/b_l;
          rhs_old = rhs[i];
          rhs[i]  = a[i]*ans[i-1] + b[i]*ans[i] + c[i]*ans[i+1];
          rhs[i]  = b_l_sq - mag + rhs[i];
          if( error > max_error )
          {
              max_error = error;
          }
       }
    }while( max_error > tol ); 
    for( i = 0; i < n_mon; i++ )
    {
       dx      = rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i];
       dy      = rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i];
       dz      = rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i];

       r1[i].x += dx;
       r1[i].y += dy;
       r1[i].z += dz;
     

     /*  v[i].x += disp1*dx/(2.0*dt);
       v[i].y += disp1*dy/(2.0*dt);
       v[i].z += disp1*dz/(2.0*dt);*/
       fc_2[i].x = 2.0*dx/(dt*dt);
       fc_2[i].y = 2.0*dy/(dt*dt);
       fc_2[i].z = 2.0*dz/(dt*dt);
    }
 
}


/* ########################################################################## */
void     constrain_positions_2filaments( vec_s  *r1,
                              vec_s  *rij1,
                              vec_s  *v1,
                              vec_s  *fc_2,
                              param_s params )
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans;

    int   n_mon, n_c, i;
    float rijx, rijy, rijz, mag;
    float rhs_old, b_l, b_l_sq, max_error, error, tol;
    float dx, dy, dz, dt;

    n_mon   = params.n_mon1;
    n_c     = n_mon - 1;
    tol     = params.tolerance;
    dt      = params.dt;
    b_l     = params.b_l;
    b_l_sq  = b_l*b_l;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_2[i].x = 0.0;
        fc_2[i].y = 0.0;
        fc_2[i].z = 0.0;
    }
    ans[0]     = 0.0;
    ans[n_mon] = 0.0;
    rij1[0].x   = 0.0;
    rij1[0].y   = 0.0;
    rij1[0].z   = 0.0;
    rij1[n_c + 1].x = 0.0;
    rij1[n_c + 1].y = 0.0;
    rij1[n_c + 1].z = 0.0;

    rijx  = r1[0].x - r1[1].x;
    rijy  = r1[0].y - r1[1].y;
    rijz  = r1[0].z - r1[1].z;
    mag   = rijx*rijx + rijy*rijy + rijz*rijz; 

    a[1]   = 0.0;
    b[1]   = rij1[1].x*rijx;
    b[1]  += rij1[1].y*rijy;
    b[1]  += rij1[1].z*rijz;
    b[1]  *= 4.0;
    c[1]   = rij1[2].x*rijx;
    c[1]  += rij1[2].y*rijy;
    c[1]  += rij1[2].z*rijz;
    c[1]  *= -2.0;
    rhs[1] = b_l_sq - mag;

    rijx    = r1[n_c - 1].x - r1[n_c].x;
    rijy    = r1[n_c - 1].y - r1[n_c].y;
    rijz    = r1[n_c - 1].z - r1[n_c].z;
    mag     = rijx*rijx + rijy*rijy + rijz*rijz; 
    a[n_c]  = rij1[n_c - 1].x*rijx;
    a[n_c] += rij1[n_c - 1].y*rijy;
    a[n_c] += rij1[n_c - 1].z*rijz;
    a[n_c] *= -2.0;
    b[n_c]  = rij1[n_c].x*rijx;
    b[n_c] += rij1[n_c].y*rijy;
    b[n_c] += rij1[n_c].z*rijz;
    b[n_c] *= 4.0;
    c[n_c]  = 0.0;
    rhs[n_c]= b_l_sq - mag;
    for( i = 2; i < n_c; i++ )
    {
        rijx   = r1[i-1].x - r1[i].x;
        rijy   = r1[i-1].y - r1[i].y;
        rijz   = r1[i-1].z - r1[i].z;
        mag    = rijx*rijx + rijy*rijy + rijz*rijz; 
        a[i]   = -2.0*(rijx*rij1[i-1].x + rijy*rij1[i-1].y + rijz*rij1[i-1].z);
        b[i]   =  4.0*(rijx*rij1[i].x + rijy*rij1[i].y + rijz*rij1[i].z);
        c[i]   = -2.0*(rijx*rij1[i+1].x + rijy*rij1[i+1].y + rijz*rij1[i+1].z);
        rhs[i] = b_l_sq - mag;
    }
           //printf(" a constrain_positions_2filaments, abans de tridag \n");

    do
    { 
           //printf(" a constrain_positions_2filaments, abans de tridag \n");

       tridag( a, b, c, rhs, ans, n_c );
       //printf(" a constrain_positions_2filaments, despres de tridag \n");
       max_error = 0.0;
       for( i = 1; i <= n_c; i++ )
       {
          rijx    = r1[i-1].x + rij1[i].x*ans[i] - rij1[i-1].x*ans[i-1] -
                   (r1[i].x + rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i]);
          rijy    = r1[i-1].y + rij1[i].y*ans[i] - rij1[i-1].y*ans[i-1] -
                   (r1[i].y + rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i]);
          rijz    = r1[i-1].z + rij1[i].z*ans[i] - rij1[i-1].z*ans[i-1] -
                   (r1[i].z + rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i]);
          mag     = rijx*rijx + rijy*rijy + rijz*rijz;
          error   = (sqrt(mag) - b_l)/b_l;
          rhs_old = rhs[i];
          rhs[i]  = a[i]*ans[i-1] + b[i]*ans[i] + c[i]*ans[i+1];
          rhs[i]  = b_l_sq - mag + rhs[i];
          if( error > max_error )
          {
              max_error = error;
          }
       }
    }while( max_error > tol ); 

               //printf(" a constrain_positions_2filaments, despres de tridag \n");

    for( i = 0; i < n_mon; i++ )
    {
       dx      = rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i];
       dy      = rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i];
       dz      = rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i];

       r1[i].x += dx;
       r1[i].y += dy;
       r1[i].z += dz;
     

       fc_2[i].x = 2.0*dx/(dt*dt);
       fc_2[i].y = 2.0*dy/(dt*dt);
       fc_2[i].z = 2.0*dz/(dt*dt);
    }
 
}


/* ########################################################################## */
void     constrain_positions_circ( vec_s  *r1,
                              vec_s  *rij1,
                              vec_s  *v1,
                              vec_s  *fc_2,
                              param_s params )
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans;

    int   n_mon, n_c, i;
    float rijx, rijy, rijz, mag;
    float rhs_old, b_l, b_l_sq, max_error, error, tol;
    float dx, dy, dz, dt, disp1;
    int indexm, indexM;
    float alpha, beta;

    n_mon   = params.n_mon;
    n_c     = n_mon - 1;
    tol     = params.tolerance;
    dt      = params.dt;
    b_l     = params.b_l;
    b_l_sq  = b_l*b_l;
    disp1    = 1.0/(1.0 + 0.5*params.gamma*params.dt);

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_2[i].x = 0.0;
        fc_2[i].y = 0.0;
        fc_2[i].z = 0.0;
    }
    //ans[0]     = 0.0;
    //ans[n_mon] = 0.0;
  
    //rij1[n_c + 1].x = 0.0;
    //rij1[n_c + 1].y = 0.0;
    //rij1[n_c + 1].z = 0.0;


// MODIFICAT FILAMENT CIRCULAR -----------
    rijx  = r1[0].x - r1[1].x;
    rijy  = r1[0].y - r1[1].y;
    rijz  = r1[0].z - r1[1].z;
    mag   = rijx*rijx + rijy*rijy + rijz*rijz; 


    beta = -2.*(rij1[n_mon].x*rijx + rij1[n_mon].y*rijy + rij1[n_mon].z*rijz);
    a[1]   = 0.;
    b[1]   = 4.0*(rij1[1].x*rijx + rij1[1].y*rijy + rij1[1].z*rijz);
    c[1]   = -2.0*(rij1[2].x*rijx + rij1[2].y*rijy +rij1[2].z*rijz);
    rhs[1] = b_l_sq - mag;
// -----------------------------------------


// MODIFICAT FILAMENT CIRCULAR -----------
    rijx    = r1[n_mon-2].x - r1[n_mon-1].x;
    rijy    = r1[n_mon-2].y - r1[n_mon-1].y;
    rijz    = r1[n_mon-2].z - r1[n_mon-1].z;
    mag     = rijx*rijx + rijy*rijy + rijz*rijz; 

    a[n_mon-1]  = rij1[n_mon-2].x*rijx;
    a[n_mon-1] += rij1[n_mon-2].y*rijy;
    a[n_mon-1] += rij1[n_mon-2].z*rijz;
    a[n_mon-1] *= -2.0;
    b[n_mon-1]  = rij1[n_mon-1].x*rijx;
    b[n_mon-1] += rij1[n_mon-1].y*rijy;
    b[n_mon-1] += rij1[n_mon-1].z*rijz;
    b[n_mon-1] *= 4.0;
    c[n_mon-1]  = rij1[n_mon].x*rijx;
    c[n_mon-1]  += rij1[n_mon].y*rijy;
    c[n_mon-1]  += rij1[n_mon].z*rijz;
    c[n_mon-1]  *= -2.0;

    rhs[n_mon-1]= b_l_sq - mag;
// -------------------------------------


// AFEGIT FILAMENT CIRCULAR ----------- 
    rijx  = r1[n_mon-1].x - r1[0].x;
    rijy  = r1[n_mon-1].y - r1[0].y;
    rijz  = r1[n_mon-1].z - r1[0].z;
    mag   = rijx*rijx + rijy*rijy + rijz*rijz; 


    a[n_mon] = -2.*(rij1[n_mon-1].x*rijx + rij1[n_mon-1].y*rijy + rij1[n_mon-1].z*rijz);
    b[n_mon] = 4.*(rij1[n_mon].x*rijx + rij1[n_mon].y*rijy + rij1[n_mon].z*rijz);
    c[n_mon] = 0.0;
    alpha = -2.*(rij1[1].x*rijx + rij1[1].y*rijy + rij1[1].z*rijz);
  
    rhs[n_mon] = b_l_sq - mag;
// -------------------------------------
    // //printf("rhs[n_mon] = %f, rij1[n_mon] = %f\n",rhs[n_mon], rij1[n_mon].x);

    for( i = 2; i < n_mon-1; i++ )
    {
        rijx   = r1[i-1].x - r1[i].x;
        rijy   = r1[i-1].y - r1[i].y;
        rijz   = r1[i-1].z - r1[i].z;
        mag    = rijx*rijx + rijy*rijy + rijz*rijz; 
        a[i]   = -2.0*(rijx*rij1[i-1].x + rijy*rij1[i-1].y + rijz*rij1[i-1].z);
        b[i]   =  4.0*(rijx*rij1[i].x + rijy*rij1[i].y + rijz*rij1[i].z);
        c[i]   = -2.0*(rijx*rij1[i+1].x + rijy*rij1[i+1].y + rijz*rij1[i+1].z);
        rhs[i] = b_l_sq - mag;
    }
    do
    {  

       cyclic( a, b, c, alpha, beta, rhs, ans, n_mon);
       //tridag( a, b, c, rhs, ans, n_mon-1);
       max_error = 0.0;
       for( i = 1; i <= n_mon; i++ )
       {

        if(i==1){
          rijx    = r1[i-1].x + rij1[i].x*ans[i] - rij1[n_mon].x*ans[n_mon] -
                   (r1[i].x + rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i]);
          rijy    = r1[i-1].y + rij1[i].y*ans[i] - rij1[n_mon].y*ans[n_mon] -
                   (r1[i].y + rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i]);
          rijz    = r1[i-1].z + rij1[i].z*ans[i] - rij1[n_mon].z*ans[n_mon] -
                   (r1[i].z + rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i]);
          mag     = rijx*rijx + rijy*rijy + rijz*rijz;
          error   = (sqrt(mag) - b_l)/b_l;
          rhs_old = rhs[i];
          rhs[i]  = a[i]*ans[n_mon] + b[i]*ans[i] + c[i]*ans[i+1];
          rhs[i]  = b_l_sq - mag + rhs[i];
       // //printf("i = %i, rhs[i] = %f, ans[i] = %f\n, mag = %f\n",i, rhs[i], ans[i], mag);

        }

        else if(i==n_mon){
          rijx    = r1[i-1].x + rij1[i].x*ans[i] - rij1[i-1].x*ans[i-1] -
                   (r1[0].x + rij1[1].x*ans[1] - rij1[i].x*ans[i]);
          rijy    = r1[i-1].y + rij1[i].y*ans[i] - rij1[i-1].y*ans[i-1] -
                   (r1[0].y + rij1[1].y*ans[1] - rij1[i].y*ans[i]);
          rijz    = r1[i-1].z + rij1[i].z*ans[i] - rij1[i-1].z*ans[i-1] -
                   (r1[0].z + rij1[1].z*ans[1] - rij1[i].z*ans[i]);
          mag     = rijx*rijx + rijy*rijy + rijz*rijz;
          error   = (sqrt(mag) - b_l)/b_l;
          rhs_old = rhs[i];
          rhs[i]  = a[i]*ans[i-1] + b[i]*ans[i] + c[i]*ans[1];
          rhs[i]  = b_l_sq - mag + rhs[i];
       // //printf("i = %i, rhs[i] = %f, ans[i] = %f\n, mag = %f\n",i, rhs[i], ans[i], mag);

        }
        else{

          rijx    = r1[i-1].x + rij1[i].x*ans[i] - rij1[i-1].x*ans[i-1] -
                   (r1[i].x + rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i]);
          rijy    = r1[i-1].y + rij1[i].y*ans[i] - rij1[i-1].y*ans[i-1] -
                   (r1[i].y + rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i]);
          rijz    = r1[i-1].z + rij1[i].z*ans[i] - rij1[i-1].z*ans[i-1] -
                   (r1[i].z + rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i]);
          mag     = rijx*rijx + rijy*rijy + rijz*rijz;
          error   = (sqrt(mag) - b_l)/b_l;
          rhs_old = rhs[i];
          rhs[i]  = a[i]*ans[i-1] + b[i]*ans[i] + c[i]*ans[i+1];
          rhs[i]  = b_l_sq - mag + rhs[i];
       // //printf("i = %i, rhs[i] = %f, ans[i] = %f\n, mag = %f\n",i, rhs[i], ans[i], mag);

        }

         if( error > max_error )
          {
              max_error = error;
          }
       }
    }while( max_error > tol ); 

    for( i = 0; i < n_mon; i++ )
    {
      if(i == 0){indexm = n_mon;}
      else{indexm = i;}

      dx      = rij1[i+1].x*ans[i+1] - rij1[indexm].x*ans[indexm];
      dy      = rij1[i+1].y*ans[i+1] - rij1[indexm].y*ans[indexm];
      dz      = rij1[i+1].z*ans[i+1] - rij1[indexm].z*ans[indexm];

       r1[i].x += dx;
       r1[i].y += dy;
       r1[i].z += dz;
     

     /*  v[i].x += disp1*dx/(2.0*dt);
       v[i].y += disp1*dy/(2.0*dt);
       v[i].z += disp1*dz/(2.0*dt);*/
       fc_2[i].x = 2.0*dx/(dt*dt);
       fc_2[i].y = 2.0*dy/(dt*dt);
       fc_2[i].z = 2.0*dz/(dt*dt);
    }
 
}
/* ########################################################################## */
void     constrain_velocities(  vec_s  *rij,
                                float  *mag,
                                vec_s  *v1,
        vec_s  *fc_1,
                                param_s params )
       
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;

    n_mon   = params.n_mon;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_1[i].x = 0.0;
        fc_1[i].y = 0.0;
        fc_1[i].z = 0.0;
    }
    for( i = 1; i <= n_c; i++ )
    {
        rhs[i]  = (v1[i-1].x - v1[i].x)*rij[i].x;
        rhs[i] += (v1[i-1].y - v1[i].y)*rij[i].y;
        rhs[i] += (v1[i-1].z - v1[i].z)*rij[i].z;
        rhs[i] /= mag[i];
        a[i]    = rij[i-1].x*rij[i].x +
                  rij[i-1].y*rij[i].y +
                  rij[i-1].z*rij[i].z;
        b[i]    = -2.0*( rij[i].x*rij[i].x +
                   rij[i].y*rij[i].y +
                   rij[i].z*rij[i].z);
        c[i]    = rij[i].x*rij[i+1].x +
                  rij[i].y*rij[i+1].y +
                  rij[i].z*rij[i+1].z;
         a[i]  /= (mag[i]*mag[i]);
         b[i]  /= (mag[i]*mag[i]);
         c[i]  /= (mag[i]*mag[i]);
     }
     tridag( a, b, c, rhs, ans, n_c );
     for( i = 1; i <= n_c; i++ )
     {

       v1[i-1].x += ans[i]*rij[i].x/mag[i];
       v1[i-1].y += ans[i]*rij[i].y/mag[i];
       v1[i-1].z += ans[i]*rij[i].z/mag[i];
       v1[i].x   -= ans[i]*rij[i].x/mag[i];
       v1[i].y   -= ans[i]*rij[i].y/mag[i];
       v1[i].z   -= ans[i]*rij[i].z/mag[i];
       fc_1[i-1].x += ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i-1].y += ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i-1].z += ans[i]*rij[i].z/(dt*mag[i]);
       fc_1[i].x   -= ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i].y   -= ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i].z   -= ans[i]*rij[i].z/(dt*mag[i]);
       
     }
     


}

/* ########################################################################## */
void     constrain_velocities_2filaments(  vec_s  *rij,
                                float  *mag,
                                vec_s  *v1,
        vec_s  *fc_1,
                                param_s params )
       
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;

    n_mon   = params.n_mon1;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_1[i].x = 0.0;
        fc_1[i].y = 0.0;
        fc_1[i].z = 0.0;
    }
    for( i = 1; i <= n_c; i++ )
    {
        rhs[i]  = (v1[i-1].x - v1[i].x)*rij[i].x;
        rhs[i] += (v1[i-1].y - v1[i].y)*rij[i].y;
        rhs[i] += (v1[i-1].z - v1[i].z)*rij[i].z;
        rhs[i] /= mag[i];
        a[i]    = rij[i-1].x*rij[i].x +
                  rij[i-1].y*rij[i].y +
                  rij[i-1].z*rij[i].z;
        b[i]    = -2.0*( rij[i].x*rij[i].x +
                   rij[i].y*rij[i].y +
                   rij[i].z*rij[i].z);
        c[i]    = rij[i].x*rij[i+1].x +
                  rij[i].y*rij[i+1].y +
                  rij[i].z*rij[i+1].z;
         a[i]  /= (mag[i]*mag[i]);
         b[i]  /= (mag[i]*mag[i]);
         c[i]  /= (mag[i]*mag[i]);
     }
     tridag( a, b, c, rhs, ans, n_c );
     for( i = 1; i <= n_c; i++ )
     {

       v1[i-1].x += ans[i]*rij[i].x/mag[i];
       v1[i-1].y += ans[i]*rij[i].y/mag[i];
       v1[i-1].z += ans[i]*rij[i].z/mag[i];
       v1[i].x   -= ans[i]*rij[i].x/mag[i];
       v1[i].y   -= ans[i]*rij[i].y/mag[i];
       v1[i].z   -= ans[i]*rij[i].z/mag[i];
       fc_1[i-1].x += ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i-1].y += ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i-1].z += ans[i]*rij[i].z/(dt*mag[i]);
       fc_1[i].x   -= ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i].y   -= ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i].z   -= ans[i]*rij[i].z/(dt*mag[i]);
       
     }
     


}


/* ########################################################################## */
void     constrain_velocities_circ(  vec_s  *rij,
                                float  *mag,
                                vec_s  *v1,
        vec_s  *fc_1,
                                param_s params )
       
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;

    int indexm;
    float alpha, beta;


    n_mon   = params.n_mon;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_1[i].x = 0.0;
        fc_1[i].y = 0.0;
        fc_1[i].z = 0.0;
    }
    for( i = 2; i <= n_mon-1; i++ )
    {

        rhs[i]  = (v1[i-1].x - v1[i].x)*rij[i].x;
        rhs[i] += (v1[i-1].y - v1[i].y)*rij[i].y;
        rhs[i] += (v1[i-1].z - v1[i].z)*rij[i].z;
        rhs[i] /= mag[i];
        a[i]    = rij[i-1].x*rij[i].x +
                  rij[i-1].y*rij[i].y +
                  rij[i-1].z*rij[i].z;
        b[i]    = -2.0*( rij[i].x*rij[i].x +
                   rij[i].y*rij[i].y +
                   rij[i].z*rij[i].z);
        c[i]    = rij[i].x*rij[i+1].x +
                  rij[i].y*rij[i+1].y +
                  rij[i].z*rij[i+1].z;
         a[i]  /= (mag[i]*mag[i]);
         b[i]  /= (mag[i]*mag[i]);
         c[i]  /= (mag[i]*mag[i]);
         // //printf("i = %i, ai = %f\n",i, a[i]);
     }

    rhs[1]  = (v1[0].x - v1[1].x)*rij[1].x;
    rhs[1] += (v1[0].y - v1[1].y)*rij[1].y;
    rhs[1] += (v1[0].z - v1[1].z)*rij[1].z;
    rhs[1] /= mag[1];
    a[1]    = 0.0;
    beta = rij[n_mon].x*rij[1].x + rij[n_mon].y*rij[1].y + rij[n_mon].z*rij[1].z;
    b[1]    = -2.0*( rij[1].x*rij[1].x + rij[1].y*rij[1].y + rij[1].z*rij[1].z);
    c[1]    = rij[1].x*rij[2].x + rij[1].y*rij[2].y + rij[1].z*rij[2].z;
     beta  /= (mag[1]*mag[1]);
     b[1]  /= (mag[1]*mag[1]);
     c[1]  /= (mag[1]*mag[1]);

    rhs[n_mon]  = (v1[n_mon-1].x - v1[0].x)*rij[n_mon].x + (v1[n_mon-1].y - v1[0].y)*rij[n_mon].y + (v1[n_mon-1].z - v1[0].z)*rij[n_mon].z;
    rhs[n_mon] /= mag[n_mon];
    // //printf("rhs[n_mon] = %f, mag[n_mon] = %f\n",rhs[n_mon], mag[n_mon]);

    a[n_mon]    = rij[n_mon-1].x*rij[n_mon].x + rij[n_mon-1].y*rij[n_mon].y + rij[n_mon-1].z*rij[n_mon].z;
    b[n_mon]    = -2.0*(rij[n_mon].x*rij[n_mon].x + rij[n_mon].y*rij[n_mon].y + rij[n_mon].z*rij[n_mon].z);
    alpha = rij[1].x*rij[n_mon].x + rij[1].y*rij[n_mon].y + rij[1].z*rij[n_mon].z;
    c[n_mon]    = 0.0;
     a[n_mon]  /= (mag[n_mon]*mag[n_mon]);
     b[n_mon]  /= (mag[n_mon]*mag[n_mon]);
     alpha  /= (mag[n_mon]*mag[n_mon]);
     // //printf("rhs[n_mon-1] = %f, mag[n_mon] = %f\n",rhs[n_mon-1], mag[n_mon]);

     cyclic( a, b, c, alpha, beta, rhs, ans, n_mon );
     //tridag(a,b,c,rhs,ans,n_mon-1);
     for( i = 0; i <= n_mon-1; i++ )
     {
      if(i==0){indexm = n_mon;}
      else{indexm = i;}
     // //printf("i = %i, rhs[i] = %f, ans[i] = %f\n",i, rhs[indexm], ans[indexm]);
      v1[i].x += ans[i+1]*rij[i+1].x/mag[i+1] -  ans[indexm]*rij[indexm].x/mag[indexm];
      v1[i].y += ans[i+1]*rij[i+1].y/mag[i+1] -  ans[indexm]*rij[indexm].y/mag[indexm];
      v1[i].z += ans[i+1]*rij[i+1].z/mag[i+1] -  ans[indexm]*rij[indexm].z/mag[indexm];

      fc_1[i].x += ans[i+1]*rij[i+1].x/(dt*mag[i+1]) -  ans[indexm]*rij[indexm].x/(dt*mag[indexm]);
      fc_1[i].y += ans[i+1]*rij[i+1].y/(dt*mag[i+1]) -  ans[indexm]*rij[indexm].y/(dt*mag[indexm]);
      fc_1[i].z += ans[i+1]*rij[i+1].z/(dt*mag[i+1]) -  ans[indexm]*rij[indexm].z/(dt*mag[indexm]);

/*       v1[i-1].x += ans[i]*rij[i].x/mag[i];
       v1[i-1].y += ans[i]*rij[i].y/mag[i];
       v1[i-1].z += ans[i]*rij[i].z/mag[i];
       v1[i].x   -= ans[i]*rij[i].x/mag[i];
       v1[i].y   -= ans[i]*rij[i].y/mag[i];
       v1[i].z   -= ans[i]*rij[i].z/mag[i];
       fc_1[i-1].x += ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i-1].y += ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i-1].z += ans[i]*rij[i].z/(dt*mag[i]);
       fc_1[i].x   -= ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i].y   -= ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i].z   -= ans[i]*rij[i].z/(dt*mag[i]);*/
       
     }
     


}

/* ########################################################################## */
void     constrain_forces(      vec_s  *rij1,
                                float  *mags1,
                                vec_s  *f,
                                vec_s  *fc_3,
                                param_s params )

/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;

    n_mon   = params.n_mon;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        b      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        c      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_3[i].x = 0.0;
        fc_3[i].y = 0.0;
        fc_3[i].z = 0.0;
    }
    for( i = 1; i <= n_c; i++ )
    {
        rhs[i]  = (f[i-1].x - f[i].x)*rij1[i].x;
        rhs[i] += (f[i-1].y - f[i].y)*rij1[i].y;
        rhs[i] += (f[i-1].z - f[i].z)*rij1[i].z;
        rhs[i] /= mags1[i];
        a[i]    = rij1[i-1].x*rij1[i].x +
                  rij1[i-1].y*rij1[i].y +
                  rij1[i-1].z*rij1[i].z;
        b[i]    = -2.0*( rij1[i].x*rij1[i].x +
                   rij1[i].y*rij1[i].y +
                   rij1[i].z*rij1[i].z);
        c[i]    = rij1[i].x*rij1[i+1].x +
                  rij1[i].y*rij1[i+1].y +
                  rij1[i].z*rij1[i+1].z;
         a[i]  /= (mags1[i]*mags1[i]);
         b[i]  /= (mags1[i]*mags1[i]);
         c[i]  /= (mags1[i]*mags1[i]);
     }
     tridag( a, b, c, rhs, ans, n_c );
     for( i = 1; i <= n_c; i++ )
     {
         fc_3[i-1].x += ans[i]*rij1[i].x/(mags1[i]);
         fc_3[i-1].y += ans[i]*rij1[i].y/(mags1[i]);
         fc_3[i-1].z += ans[i]*rij1[i].z/(mags1[i]);
         fc_3[i].x   -= ans[i]*rij1[i].x/(mags1[i]);
         fc_3[i].y   -= ans[i]*rij1[i].y/(mags1[i]);
         fc_3[i].z   -= ans[i]*rij1[i].z/(mags1[i]);

     }
}

/* ########################################################################## */
void     constrain_forces_2filaments(      vec_s  *rij1,
                                float  *mags1,
                                vec_s  *f,
                                vec_s  *fc_3,
                                param_s params )

/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;

    n_mon   = params.n_mon1;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        b      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        c      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_3[i].x = 0.0;
        fc_3[i].y = 0.0;
        fc_3[i].z = 0.0;
    }
    for( i = 1; i <= n_c; i++ )
    {
        rhs[i]  = (f[i-1].x - f[i].x)*rij1[i].x;
        rhs[i] += (f[i-1].y - f[i].y)*rij1[i].y;
        rhs[i] += (f[i-1].z - f[i].z)*rij1[i].z;
        rhs[i] /= mags1[i];
        a[i]    = rij1[i-1].x*rij1[i].x +
                  rij1[i-1].y*rij1[i].y +
                  rij1[i-1].z*rij1[i].z;
        b[i]    = -2.0*( rij1[i].x*rij1[i].x +
                   rij1[i].y*rij1[i].y +
                   rij1[i].z*rij1[i].z);
        c[i]    = rij1[i].x*rij1[i+1].x +
                  rij1[i].y*rij1[i+1].y +
                  rij1[i].z*rij1[i+1].z;
         a[i]  /= (mags1[i]*mags1[i]);
         b[i]  /= (mags1[i]*mags1[i]);
         c[i]  /= (mags1[i]*mags1[i]);
     }
     tridag( a, b, c, rhs, ans, n_c );
     for( i = 1; i <= n_c; i++ )
     {
         fc_3[i-1].x += ans[i]*rij1[i].x/(mags1[i]);
         fc_3[i-1].y += ans[i]*rij1[i].y/(mags1[i]);
         fc_3[i-1].z += ans[i]*rij1[i].z/(mags1[i]);
         fc_3[i].x   -= ans[i]*rij1[i].x/(mags1[i]);
         fc_3[i].y   -= ans[i]*rij1[i].y/(mags1[i]);
         fc_3[i].z   -= ans[i]*rij1[i].z/(mags1[i]);

     }
}


/* ########################################################################## */
void     constrain_forces_circ(      vec_s  *rij1,
                                float  *mag,
                                vec_s  *f,
                                vec_s  *fc_3,
                                param_s params )

/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;
    int indexm;
    float alpha, beta;

    n_mon   = params.n_mon;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        b      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        c      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_3[i].x = 0.0;
        fc_3[i].y = 0.0;
        fc_3[i].z = 0.0;
    }

    for( i = 2; i <= n_mon-1; i++ )
    {
        rhs[i]  = (f[i-1].x - f[i].x)*rij1[i].x;
        rhs[i] += (f[i-1].y - f[i].y)*rij1[i].y;
        rhs[i] += (f[i-1].z - f[i].z)*rij1[i].z;
        rhs[i] /= mag[i];
        a[i]    = rij1[i-1].x*rij1[i].x +
                  rij1[i-1].y*rij1[i].y +
                  rij1[i-1].z*rij1[i].z;
        b[i]    = -2.0*( rij1[i].x*rij1[i].x +
                   rij1[i].y*rij1[i].y +
                   rij1[i].z*rij1[i].z);
        c[i]    = rij1[i].x*rij1[i+1].x +
                  rij1[i].y*rij1[i+1].y +
                  rij1[i].z*rij1[i+1].z;
         a[i]  /= (mag[i]*mag[i]);
         b[i]  /= (mag[i]*mag[i]);
         c[i]  /= (mag[i]*mag[i]);
     }

    rhs[1]  = (f[0].x - f[1].x)*rij1[1].x;
    rhs[1] += (f[0].y - f[1].y)*rij1[1].y;
    rhs[1] += (f[0].z - f[1].z)*rij1[1].z;
    rhs[1] /= mag[1];
    a[1]    = 0.0;
    beta = rij1[n_mon].x*rij1[1].x + rij1[n_mon].y*rij1[1].y + rij1[n_mon].z*rij1[1].z;
    b[1]    = -2.0*( rij1[1].x*rij1[1].x + rij1[1].y*rij1[1].y + rij1[1].z*rij1[1].z);
    c[1]    = rij1[1].x*rij1[2].x + rij1[1].y*rij1[2].y + rij1[1].z*rij1[2].z;
     beta  /= (mag[1]*mag[1]);
     b[1]  /= (mag[1]*mag[1]);
     c[1]  /= (mag[1]*mag[1]);

    rhs[n_mon]  = (f[n_mon-1].x - f[0].x)*rij1[n_mon].x + (f[n_mon-1].y - f[0].y)*rij1[n_mon].y + (f[n_mon-1].z - f[0].z)*rij1[n_mon].z;
    rhs[n_mon] /= mag[n_mon];
    a[n_mon]    = rij1[n_mon-1].x*rij1[n_mon].x + rij1[n_mon-1].y*rij1[n_mon].y + rij1[n_mon-1].z*rij1[n_mon].z;
    b[n_mon]    = -2.0*(rij1[n_mon].x*rij1[n_mon].x + rij1[n_mon].y*rij1[n_mon].y + rij1[n_mon].z*rij1[n_mon].z);
    alpha = rij1[1].x*rij1[n_mon].x + rij1[1].y*rij1[n_mon].y + rij1[1].z*rij1[n_mon].z;
    c[n_mon]    = 0.0;
     a[n_mon]  /= (mag[n_mon]*mag[n_mon]);
     b[n_mon]  /= (mag[n_mon]*mag[n_mon]);
     alpha  /= (mag[n_mon]*mag[n_mon]);


     cyclic( a, b, c, alpha, beta, rhs, ans, n_mon );
     for( i = 0; i <= n_mon-1; i++ )
     {
      if(i==0){indexm = n_mon;}
      else{indexm = i;}

      fc_3[i].x += ans[i+1]*rij1[i+1].x/mag[i+1] -  ans[indexm]*rij1[indexm].x/mag[indexm];
      fc_3[i].y += ans[i+1]*rij1[i+1].y/mag[i+1] -  ans[indexm]*rij1[indexm].y/mag[indexm];
      fc_3[i].z += ans[i+1]*rij1[i+1].z/mag[i+1] -  ans[indexm]*rij1[indexm].z/mag[indexm];
      
      }

}

/* ########################################################################## */
float    verlet_pt1(  vec_s  *r,
                      vec_s  *v,
                      vec_s  *fa,
                      param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq;
  static float hydr_radius, gamma, massa;
  static float disp1, disp2;
  static float bond;
  static vec_s *v_h = NULL;
  static int  n_mon, n_fil;

  float ddx, ddy, ddz;
  int  i;
  vec_s f_i, v_norm, v_hydr;
  

    if( init_flag )
    {

        dt             = parameters.dt;
  bond           = parameters.b_l;
  hydr_radius    = parameters.hydr_radius;
  gamma          = parameters.gamma /*  * (hydr_radius*bond) */;
  massa          = parameters.mass;
  disp1          = 1.0/(1.0 + 0.5*gamma*dt);
  n_fil          = parameters.n_fil;
  n_mon          = parameters.n_mon;

        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
        init_flag  = FALSE;
    }

    if(v_h == NULL){
      v_h = (vec_s *)malloc(n_mon*n_fil*sizeof( vec_s ));     
    }

    if (parameters.oseen ==1){
      oseen(r, fa, v_h, parameters);
    }
    if (parameters.blake ==1){
      blake(r, fa, v_h, parameters);
    }

   
    for( i = 0; i < n_mon*n_fil; i++ )
    {

        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;

  
        v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
    
        v_hydr.x   = v_h[i].x; 
        v_hydr.y   = v_h[i].y;
        v_hydr.z   = v_h[i].z;

        ddx     = v_norm.x*dt + (f_i.x + gamma*(v_hydr.x-v_norm.x))*half_dt_sq;
        ddy     = v_norm.y*dt + (f_i.y + gamma*(v_hydr.y-v_norm.y))*half_dt_sq;
        ddz     = v_norm.z*dt + (f_i.z + gamma*(v_hydr.z-v_norm.z))*half_dt_sq;


        r[i].x      += ddx;
        r[i].y      += ddy;
        r[i].z      += ddz;
   
  

     v_norm.x = disp1*(v_norm.x + (f_i.x + gamma*(v_hydr.x-v_norm.x))*half_dt);
     v_norm.y = disp1*(v_norm.y + (f_i.y + gamma*(v_hydr.y-v_norm.y))*half_dt);
     v_norm.z = disp1*(v_norm.z + (f_i.z + gamma*(v_hydr.z-v_norm.z))*half_dt);


       

       v[i].x        = v_norm.x;
       v[i].y        = v_norm.y;
       v[i].z        = v_norm.z;



    }
    
   
    return( parameters.mass/2.0);
} 

/* ########################################################################## */
float    verlet_pt1_2filaments(  vec_s  *r,
                      vec_s  *v,
                      vec_s  *fa,
                      param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq;
  static float hydr_radius, gamma;
  static float disp1;
  static float bond;
  static vec_s *v_h = NULL;
  static int  n_mon1, n_mon2, n_fil;

  float ddx, ddy, ddz, f_filx;
  int  i;
  vec_s f_i, v_norm, v_hydr;
  

    if( init_flag )
    {

        dt             = parameters.dt;
  bond           = parameters.b_l;
  hydr_radius    = parameters.hydr_radius;
  n_fil          = parameters.n_fil;
  n_mon1          = parameters.n_mon1;
  n_mon2          = parameters.n_mon2;

        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
        init_flag  = FALSE;
    }

    if(v_h == NULL){
      v_h = (vec_s *)malloc((n_mon1+n_mon2)*sizeof( vec_s ));     
    }

    if (parameters.oseen ==1){
      oseen_2filaments(r, fa, v_h, parameters);
    }
    if (parameters.blake ==1){
      blake_2filaments(r, fa, v_h, parameters);
    }

   //printf("0: r[n_mon1].x = %le\n", r[n_mon1].x );
   //printf("0: r[0].x = %le\n", r[0].x );
    f_filx = 0.0;
    for( i = 0; i < (n_mon1 + n_mon2); i++ )
    {
      if(i<n_mon1){gamma = parameters.gamma1;
        f_filx += fa[i].x*parameters.mass1;
      }
      else{gamma = parameters.gamma2;}

        disp1         = 1.0/(1.0 + 0.5*gamma*dt);

        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;

  
        v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
    
        v_hydr.x   = v_h[i].x; 
        v_hydr.y   = v_h[i].y;
        v_hydr.z   = v_h[i].z;

        ddx     = v_norm.x*dt + (f_i.x + gamma*(v_hydr.x-v_norm.x))*half_dt_sq;
        ddy     = v_norm.y*dt + (f_i.y + gamma*(v_hydr.y-v_norm.y))*half_dt_sq;
        ddz     = v_norm.z*dt + (f_i.z + gamma*(v_hydr.z-v_norm.z))*half_dt_sq;


        r[i].x      += ddx;
        r[i].y      += ddy;
        r[i].z      += ddz;

        v_norm.x = disp1*(v_norm.x + (f_i.x + gamma*(v_hydr.x-v_norm.x))*half_dt);
        v_norm.y = disp1*(v_norm.y + (f_i.y + gamma*(v_hydr.y-v_norm.y))*half_dt);
        v_norm.z = disp1*(v_norm.z + (f_i.z + gamma*(v_hydr.z-v_norm.z))*half_dt);

        v[i].x        = v_norm.x;
        v[i].y        = v_norm.y;
        v[i].z        = v_norm.z;
                    //printf("%le\t%le\n", v[i].x, r[i].x);
//if(i==n_mon1){printf("f_filx = %le, fx = %le, fv = %le,  f - fv = %le,  vx = %le,  vhydrx = %le\n", f_filx, f_i.x*parameters.mass2, gamma*(v_hydr.x-v_norm.x), f_i.x + gamma*(v_hydr.x-v_norm.x), v_norm.x, v_hydr.x);}

    }
  
   //printf("1: r[n_mon1].x = %le\n", r[n_mon1].x );
   //printf("1: r[0].x = %le\n", r[0].x );

    return( parameters.mass/2.0);
} 


/* ########################################################################## */
float    verlet_pt1_gamma_var(  vec_s  *r,
                      vec_s  *v,
                      vec_s  *fa,
                      param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq;
  static float hydr_radius, gamma, massa;
  static float disp1, disp2;
  static float bond;
  static vec_s *v_h = NULL;
  static int  n_mon, n_fil;

  float ddx, ddy, ddz, gamma1, gamma2, gam;
  int  i;
  vec_s f_i, v_norm, v_hydr;
  

    if( init_flag )
    {

  dt             = parameters.dt;
  bond           = parameters.b_l;
  hydr_radius    = parameters.hydr_radius;
  gamma          = parameters.gamma /*  * (hydr_radius*bond) */;
  massa          = parameters.mass;
  disp1          = 1.0/(1.0 + 0.5*gamma*dt);
  n_fil          = parameters.n_fil;
  n_mon          = parameters.n_mon;

  gamma1         = gamma*10.; // gamma paramagnet
  gamma2         = gamma*100.; // gamma ferromagnet


        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
        init_flag  = FALSE;
    }

    if(v_h == NULL){
      v_h = (vec_s *)malloc(n_mon*n_fil*sizeof( vec_s ));     
    }

    if (parameters.oseen ==1){
      oseen(r, fa, v_h, parameters);
    }
    if (parameters.blake ==1){
      blake(r, fa, v_h, parameters);
    }

   
    for( i = 0; i < n_mon*n_fil; i++ )
    {
        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;
  
        v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
 
        v_hydr.x   = v_h[i].x; 
        v_hydr.y   = v_h[i].y;
        v_hydr.z   = v_h[i].z;

  if(i==0){gam = gamma1;}
  else if(i==n_mon-1){gam = gamma2;}
  else{gam = gamma;}
  disp1    = 1.0/(1.0 + 0.5*gam*dt);

  ddx     = v_norm.x*dt + (f_i.x + gam*(v_hydr.x-v_norm.x))*half_dt_sq;
  ddy     = v_norm.y*dt + (f_i.y + gam*(v_hydr.y-v_norm.y))*half_dt_sq;
  ddz     = v_norm.z*dt + (f_i.z + gam*(v_hydr.z-v_norm.z))*half_dt_sq;
  r[i].x      += ddx;
  r[i].y      += ddy;
  r[i].z      += ddz;

   v_norm.x = disp1*(v_norm.x + (f_i.x + gam*(v_hydr.x-v_norm.x))*half_dt);
   v_norm.y = disp1*(v_norm.y + (f_i.y + gam*(v_hydr.y-v_norm.y))*half_dt);
   v_norm.z = disp1*(v_norm.z + (f_i.z + gam*(v_hydr.z-v_norm.z))*half_dt);

   v[i].x        = v_norm.x;
   v[i].y        = v_norm.y;
   v[i].z        = v_norm.z;

    }
    
   
    return( parameters.mass/2.0);
} 




/* ########################################################################## */
float    verlet_pt1_filament(  vec_s  *r,
                      vec_s  *v,
                      vec_s  *fa,
                      param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq;
  static float hydr_radius, gamma, massa;
  static float disp1, disp2;
  static float bond;
  static vec_s *v_h = NULL;
  static int  n_mon, n_fil;

  float ddx, ddy, ddz, gammapara, gammaferro, gam;
  int  i;
  vec_s f_i, v_norm, v_hydr;
  

    if( init_flag )
    {

  dt             = parameters.dt;
  bond           = parameters.b_l;
  hydr_radius    = parameters.hydr_radius;
  gamma          = parameters.gamma /*  * (hydr_radius*bond) */;
//  gammapara         = parameters.gamma1; // gamma paramagnet
//  gammaferro         = parameters.gamma2; // gamma ferromagnet
  gammapara         = parameters.gamma1; // gamma paramagnet
  gammaferro        = parameters.gamma2; // gamma ferromagnet
  massa          = parameters.mass;
  disp1          = 1.0/(1.0 + 0.5*gamma*dt);
  n_fil          = parameters.n_fil;
  n_mon          = parameters.n_mon;


        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
        init_flag  = FALSE;
    }

    if(v_h == NULL){
      v_h = (vec_s *)malloc(n_mon*n_fil*sizeof( vec_s ));     
    }

    if (parameters.oseen ==1){
      oseen(r, fa, v_h, parameters);
    }
    if (parameters.blake ==1){
      blake(r, fa, v_h, parameters);
    }

   
    for( i = 0; i < n_mon*n_fil; i++ )
    {
        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;
  
        v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
 
        v_hydr.x   = v_h[i].x; 
        v_hydr.y   = v_h[i].y;
        v_hydr.z   = v_h[i].z;

  if(i==0){gam = gammapara;}
  else if(i==n_mon-1){gam = gammaferro;}
  else{gam = gamma;}
  disp1    = 1.0/(1.0 + 0.5*gam*dt);

  ddx     = v_norm.x*dt + (f_i.x + gam*(v_hydr.x-v_norm.x))*half_dt_sq;
  ddy     = v_norm.y*dt + (f_i.y + gam*(v_hydr.y-v_norm.y))*half_dt_sq;
  ddz     = v_norm.z*dt + (f_i.z + gam*(v_hydr.z-v_norm.z))*half_dt_sq;
  r[i].x      += ddx;
  r[i].y      += ddy;
  r[i].z      += ddz;

   v_norm.x = disp1*(v_norm.x + (f_i.x + gam*(v_hydr.x-v_norm.x))*half_dt);
   v_norm.y = disp1*(v_norm.y + (f_i.y + gam*(v_hydr.y-v_norm.y))*half_dt);
   v_norm.z = disp1*(v_norm.z + (f_i.z + gam*(v_hydr.z-v_norm.z))*half_dt);

   v[i].x        = v_norm.x;
   v[i].y        = v_norm.y;
   v[i].z        = v_norm.z;

    }
    
   
    return( parameters.mass/2.0);
} 

/* ########################################################################## */
void      verlet_pt2(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *fa,
                       param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq, disp1;
  static float gamma, massa;
  static float bond, hydr_radius;
  static vec_s *v_hyd = NULL;
  static int n_mon, n_fil;
  int  i;

  vec_s f_i, v_norm, v_hydro;

    if( init_flag )
    {
        dt         = parameters.dt;
	bond       = parameters.b_l;
	hydr_radius = parameters.hydr_radius;
	gamma      = parameters.gamma;
        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
	massa      = parameters.mass;
	disp1      = 1.0/(1.0 + 0.5*gamma*dt);
	n_fil      = parameters.n_fil;
	n_mon      = parameters.n_mon;

        init_flag  = FALSE;
    }
    if(v_hyd == NULL)
      v_hyd = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));

    if (parameters.oseen ==1){
      oseen(r, fa, v_hyd, parameters);
    }
    if (parameters.blake ==1){
      blake(r, fa, v_hyd, parameters);
    }

    
    
    for( i = 0; i < n_fil*n_mon; i++ )
    {

        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;
   

	      v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
	
	v_hydro.x   = v_hyd[i].x;
 	v_hydro.y   = v_hyd[i].y;
	v_hydro.z   = v_hyd[i].z;


     v_norm.x = disp1* (f_i.x + gamma*(v_hydro.x))*half_dt;
     v_norm.y = disp1* (f_i.y + gamma*(v_hydro.y))*half_dt;
     v_norm.z = disp1* (f_i.z + gamma*(v_hydro.z))*half_dt;
       
       v[i].x      += v_norm.x;
       v[i].y      += v_norm.y;
       v[i].z      += v_norm.z;


    }
}


/* ########################################################################## */
void      verlet_pt2_2filaments(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *fa,
                       param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq, disp1;
  static float gamma;
  static float bond, hydr_radius;
  static vec_s *v_hyd = NULL;
  static int n_mon1, n_mon2, n_fil;
  int  i;

  vec_s f_i, v_norm, v_hydro;

    if( init_flag )
    {
        dt         = parameters.dt;
        bond       = parameters.b_l;
        hydr_radius = parameters.hydr_radius;
        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
        n_fil      = parameters.n_fil;
        n_mon1      = parameters.n_mon1;
        n_mon2      = parameters.n_mon2;

        init_flag  = FALSE;
    }
    if(v_hyd == NULL)
      v_hyd = (vec_s *)malloc( (n_mon1+n_mon2)*sizeof( vec_s ));

    if (parameters.oseen ==1){
      oseen_2filaments(r, fa, v_hyd, parameters);
    }
    if (parameters.blake ==1){
      blake_2filaments(r, fa, v_hyd, parameters);
    }
  
    for( i = 0; i < (n_mon1+n_mon2); i++ )
    {
      if(i<n_mon1){gamma = parameters.gamma1;}
      else{gamma = parameters.gamma2;}

        disp1         = 1.0/(1.0 + 0.5*gamma*dt);
        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;
   

        v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
  
        v_hydro.x   = v_hyd[i].x;
        v_hydro.y   = v_hyd[i].y;
        v_hydro.z   = v_hyd[i].z;


        v_norm.x = disp1* (f_i.x + gamma*(v_hydro.x))*half_dt;
        v_norm.y = disp1* (f_i.y + gamma*(v_hydro.y))*half_dt;
        v_norm.z = disp1* (f_i.z + gamma*(v_hydro.z))*half_dt;
       
        v[i].x      += v_norm.x;
        v[i].y      += v_norm.y;
        v[i].z      += v_norm.z;
    }
}

/* ########################################################################## */
void      verlet_pt2_gamma_var(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *fa,
                       param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq, disp1;
  static float gamma, massa;
  static float bond, hydr_radius;
  static vec_s *v_hyd = NULL;
  static int n_mon, n_fil;
  int  i;
  float gam, gamma1, gamma2;

  vec_s f_i, v_norm, v_hydro;

    if( init_flag )
    {
        dt         = parameters.dt;
  bond       = parameters.b_l;
  hydr_radius = parameters.hydr_radius;
  gamma      = parameters.gamma;
        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
  massa      = parameters.mass;

  n_fil      = parameters.n_fil;
  n_mon      = parameters.n_mon;
  gamma1 = 100.0*gamma;
  gamma2 = 10.0*gamma;
        init_flag  = FALSE;
    }
    if(v_hyd == NULL)
      v_hyd = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));

    if (parameters.oseen ==1){
      oseen(r, fa, v_hyd, parameters);
    }
    if (parameters.blake ==1){
      blake(r, fa, v_hyd, parameters);
    }

    
    
    for( i = 0; i < n_fil*n_mon; i++ )
    {

        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;
   

  v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
  
  v_hydro.x   = v_hyd[i].x;
  v_hydro.y   = v_hyd[i].y;
  v_hydro.z   = v_hyd[i].z;


  if(i==0){gam = gamma1;}
  else if(i==n_mon-1){gam = gamma2;}
  else{gam = gamma;}

     disp1    = 1.0/(1.0 + 0.5*gam*dt);
     v_norm.x = disp1* (f_i.x + gam*(v_hydro.x))*half_dt;
     v_norm.y = disp1* (f_i.y + gam*(v_hydro.y))*half_dt;
     v_norm.z = disp1* (f_i.z + gam*(v_hydro.z))*half_dt;
       
       v[i].x      += v_norm.x;
       v[i].y      += v_norm.y;
       v[i].z      += v_norm.z;


    }
}

/* ########################################################################## */
void      verlet_pt2_filament(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *fa,
                       param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq, disp1;
  static float gamma, massa;
  static float bond, hydr_radius;
  static vec_s *v_hyd = NULL;
  static int n_mon, n_fil;
  int  i;
  float gam, gammapara, gammaferro;

  vec_s f_i, v_norm, v_hydro;

    if( init_flag )
    {
      dt         = parameters.dt;
      bond       = parameters.b_l;
      hydr_radius = parameters.hydr_radius;
      gamma      = parameters.gamma;
//      gamma1         = parameters.gamma1; // gamma paramagnet
//      gamma2         = parameters.gamma2; // gamma ferromagnet
      gammapara         = parameters.gamma; // gamma paramagnet
      gammaferro         = parameters.gamma; // gamma ferromagnet
      half_dt    = dt/2.0;
      half_dt_sq = dt*dt/2.0;
      massa      = parameters.mass;

      n_fil      = parameters.n_fil;
      n_mon      = parameters.n_mon;
      init_flag  = FALSE;
    }
    if(v_hyd == NULL)
      v_hyd = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));

    if (parameters.oseen ==1){
      oseen(r, fa, v_hyd, parameters);
    }
    if (parameters.blake ==1){
      blake(r, fa, v_hyd, parameters);
    }

    
    
    for( i = 0; i < n_fil*n_mon; i++ )
    {

        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;
   

  v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
  
  v_hydro.x   = v_hyd[i].x;
  v_hydro.y   = v_hyd[i].y;
  v_hydro.z   = v_hyd[i].z;


  if(i==0){gam = gammapara;}
  else if(i==n_mon-1){gam = gammaferro;}
  else{gam = gamma;}

     disp1    = 1.0/(1.0 + 0.5*gam*dt);
     v_norm.x = disp1* (f_i.x + gam*(v_hydro.x))*half_dt;
     v_norm.y = disp1* (f_i.y + gam*(v_hydro.y))*half_dt;
     v_norm.z = disp1* (f_i.z + gam*(v_hydro.z))*half_dt;
       
       v[i].x      += v_norm.x;
       v[i].y      += v_norm.y;
       v[i].z      += v_norm.z;


    }
}

/*  ************************************************************************ */

void     oseen( vec_s *r,
		vec_s *fa,
		vec_s *v_h,
		param_s parameters)
/*  ************************************************************************ */
{
  static float hydr_radius, gamma,bond, massa, eta;
  static char  init_flag = TRUE;
  static int n_mon, n_fil;
  
  float rijx, rijy,rijz;
  float mag_ij;
  int  i,j,k; 
  vec_s r_i, r_j, f_i, f_j, v_hydr;

  if( init_flag )
    {
      eta            = parameters.viscosity;
      gamma          = parameters.gamma*parameters.mass;
      hydr_radius    = parameters.hydr_radius;
      bond           = parameters.b_l;
      massa          = parameters.mass;
      n_fil          = parameters.n_fil;
      n_mon          = parameters.n_mon;        
      init_flag  = FALSE;
    }



  /*    //printf("hydrodynamic radius is %f", hydr_radius); */
  
  for( i = 0; i < n_mon*n_fil; i++ )
    {

      v_hydr.x = 0;
      v_hydr.y = 0;
      v_hydr.z = 0;
      
      
      r_i.x      = r[i].x;
      r_i.y      = r[i].y;
      r_i.z      = r[i].z; 

      f_i.x      = massa * fa[i].x;
      f_i.y      = massa * fa[i].y;
      f_i.z      = massa * fa[i].z;
      
      for( j = 0; j <  n_mon*n_fil; j++)
	{
	  
	  f_j.x      = massa * fa[j].x;
	  f_j.y      = massa * fa[j].y;
	  f_j.z      = massa * fa[j].z;
	  
	  r_j.x      = r[j].x;
	  r_j.y      = r[j].y;
	  r_j.z      = r[j].z;
	
	  rijx     = r[i].x - r[j].x;
	  rijy     = r[i].y - r[j].y;
	  rijz     = r[i].z - r[j].z;

	  mag_ij  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
	  
	  rijx    /= mag_ij;
	  rijy    /= mag_ij;
	  rijz    /= mag_ij;
	  
//	  mag_ij  /= bond;
	  
	  if(j != i)
	    {
	      v_hydr.x += ((1+ rijx * rijx)* f_j.x +
			   (rijx * rijy)   * f_j.y +
			   (rijx * rijz)   * f_j.z  ) / mag_ij;
	      
	      v_hydr.y += ((rijx* rijy)    * f_j.x +
			   (1+ rijy * rijy)* f_j.y +
			   (rijy * rijz)   * f_j.z  ) / mag_ij;
		
	      v_hydr.z += ((rijz * rijx)   * f_j.x +
			   (rijz * rijy)   * f_j.y +
			   (1+ rijz * rijz)* f_j.z  ) / mag_ij;
	    }
	  
	  
	}

      v_h[i].x = ( v_hydr.x ) *  (1./(8.*PI*eta));
      v_h[i].y = ( v_hydr.y ) *  (1./(8.*PI*eta));
      v_h[i].z = ( v_hydr.z ) *  (1./(8.*PI*eta));
  //    v_h[i].x = (0.75 * v_hydr.x ) *  (hydr_radius/ gamma);
  //    v_h[i].y = (0.75 * v_hydr.y ) *  (hydr_radius/ gamma);
  //    v_h[i].z = (0.75 * v_hydr.z ) *  (hydr_radius/ gamma);
	
    }
}


    static  float *dx, *dy;
    static  vec_s *f_old;
    static float w;

/*  ************************************************************************ */

void     oseen_2filaments( vec_s *r,
    vec_s *fa,
    vec_s *v_h,
    param_s parameters)
/*  ************************************************************************ */
{
  static float hydr_radius, eta, gamma, bond, massa, massa1, massa2;
  static char  init_flag = TRUE;
  static int n_mon1, n_mon2, n_fil;
  
  float rijx, rijy,rijz;
  float mag_ij;
  int  i,j,k; 
  vec_s r_i, r_j, f_i, f_j, v_hydr;

  if( init_flag )
    {
      eta            = parameters.viscosity;
      gamma          = parameters.gamma1*parameters.mass1;
      massa1          = parameters.mass1;
      massa2          = parameters.mass2;
      n_fil          = parameters.n_fil;
      n_mon1          = parameters.n_mon1;        
      n_mon2          = parameters.n_mon2;        
      init_flag  = FALSE;
    }

  //printf("n_mon1 = %i, n_mon2 = %i \n",n_mon1, n_mon2);
  
  for( i = 0; i < (n_mon1 + n_mon2); i++ )
    {
      v_hydr.x = 0;
      v_hydr.y = 0;
      v_hydr.z = 0;
      
      r_i.x      = r[i].x;
      r_i.y      = r[i].y;
      r_i.z      = r[i].z; 
      

      for( j = 0; j < (n_mon1 + n_mon2); j++)
  {
    
      if(j<n_mon1){
        massa = massa1;
      }
      else{
        massa = massa2;        
      }


    f_j.x      = massa * fa[j].x;
    f_j.y      = massa * fa[j].y;
    f_j.z      = massa * fa[j].z;
    
    r_j.x      = r[j].x;
    r_j.y      = r[j].y;
    r_j.z      = r[j].z;
  
    rijx     = r[i].x - r[j].x;
    rijy     = r[i].y - r[j].y;
    rijz     = r[i].z - r[j].z;

    mag_ij  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
    
    rijx    /= mag_ij;
    rijy    /= mag_ij;
    rijz    /= mag_ij;
    
    
    if(j != i)
      {
        v_hydr.x += ((1+ rijx * rijx)* f_j.x +
         (rijx * rijy)   * f_j.y +
         (rijx * rijz)   * f_j.z  ) / mag_ij;
        
        v_hydr.y += ((rijx* rijy)    * f_j.x +
         (1+ rijy * rijy)* f_j.y +
         (rijy * rijz)   * f_j.z  ) / mag_ij;
    
        v_hydr.z += ((rijz * rijx)   * f_j.x +
         (rijz * rijy)   * f_j.y +
         (1+ rijz * rijz)* f_j.z  ) / mag_ij;
        if(i>n_mon1-2){
                  //printf("oseen %i\t%i\t%le\t%le\t%le\t%le\n",i,j,  rijx, f_j.x, v_hydr.x, r_j.x);
                }
      } 
  

  }

      v_h[i].x = ( v_hydr.x ) *  (1./(8.*PI*eta));
      v_h[i].y = ( v_hydr.y ) *  (1./(8.*PI*eta));
      v_h[i].z = ( v_hydr.z ) *  (1./(8.*PI*eta));

//      v_h[i].x = (0.75 * v_hydr.x ) *  (hydr_radius/ gamma);
//      v_h[i].y = (0.75 * v_hydr.y ) *  (hydr_radius/ gamma);
//      v_h[i].z = (0.75 * v_hydr.z ) *  (hydr_radius/ gamma);
                  //printf("%i\t%i\t%le\t%le\t%le\n",i,j, rijx, v_h[i].x, r_i.x);

    }
}

/*  ************************************************************************ */

void     blake( vec_s *r,
		vec_s *fa,
		vec_s *v_h,
		param_s parameters)
/*  ************************************************************************ */
{
  static float hydr_radius, gamma, bond, massa, eta;
  static char  init_flag = TRUE;
  static int n_mon, n_fil;
  
  float rijx, rijy,rijz,rijxbar, rijybar,rijzbar;
  float mag_ij, mag_ijbar, mag_ijbar2, mag_ijbar3, mag_ijbar5;
  int  i,j,k; 
  vec_s r_i, r_j, f_i, f_j, v_hydr;

  if( init_flag )
    {
      eta            = parameters.viscosity;
      gamma          = parameters.gamma*parameters.mass;
      hydr_radius    = parameters.hydr_radius;
      bond           = parameters.b_l;
      massa          = parameters.mass;
      n_fil          = parameters.n_fil;
      n_mon          = parameters.n_mon;        
      init_flag  = FALSE;
    }



  /*    //printf("hydrodynamic radius is %f", hydr_radius); */
  
  for( i = 0; i < n_mon*n_fil; i++ )
    {

      v_hydr.x = 0;
      v_hydr.y = 0;
      v_hydr.z = 0;
      
      
      r_i.x      = r[i].x;
      r_i.y      = r[i].y;
      r_i.z      = r[i].z; 

      f_i.x      = massa * fa[i].x;
      f_i.y      = massa * fa[i].y;
      f_i.z      = massa * fa[i].z;
      
      for( j = 0; j <  n_mon*n_fil; j++)
	{
	  
	  f_j.x      = massa * fa[j].x;
	  f_j.y      = massa * fa[j].y;
	  f_j.z      = massa * fa[j].z;
	  
	  r_j.x      = r[j].x;
	  r_j.y      = r[j].y;
	  r_j.z      = r[j].z;
	
	  rijx     = r[i].x - r[j].x;
	  rijy     = r[i].y - r[j].y;
	  rijz     = r[i].z - r[j].z;
    rijxbar     = r[i].x - r[j].x;
    rijybar     = r[i].y - r[j].y;
    rijzbar     = r[i].z + r[j].z;

	  mag_ij  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
	  mag_ijbar2  = rijx*rijx + rijy*rijy + rijzbar*rijzbar;
	  mag_ijbar  = sqrt(mag_ijbar2);
	  mag_ijbar3 = mag_ijbar*mag_ijbar2;
	  mag_ijbar5 = mag_ijbar3*mag_ijbar2;
	  
	  rijx    /= mag_ij;
	  rijy    /= mag_ij;
	  rijz    /= mag_ij;
    rijxbar  /= mag_ijbar;
    rijybar  /= mag_ijbar;
    rijzbar  /= mag_ijbar;	  
//	  mag_ij  /= bond;
//	  mag_ijbar  /= bond;
	  
	  
	  
	  if(j != i)
	    {
	    // OSEEN:
	      v_hydr.x +=  ((1+ rijx * rijx)* f_j.x +
			   (rijx * rijy)   * f_j.y +
			   (rijx * rijz)   * f_j.z  ) / mag_ij;
	      
	      v_hydr.y +=  ((rijx* rijy)    * f_j.x +
			   (1+ rijy * rijy)* f_j.y +
			   (rijy * rijz)   * f_j.z  ) / mag_ij;
		
	      v_hydr.z +=  ((rijz * rijx)   * f_j.x +
			   (rijz * rijy)   * f_j.y +
			   (1+ rijz * rijz)* f_j.z  ) / mag_ij;
			   
	    // OSEEN_BAR:
	      v_hydr.x -= ((1+ rijxbar * rijxbar)* f_j.x +
			   (rijxbar * rijybar)   * f_j.y +
			   (rijxbar * rijzbar)   * f_j.z  ) / mag_ijbar;
	      
	      v_hydr.y -= ((rijxbar* rijybar)    * f_j.x +
			   (1+ rijybar * rijybar)* f_j.y +
			   (rijybar * rijzbar)   * f_j.z  ) / mag_ijbar;
		
	      v_hydr.z -= ((rijzbar * rijxbar)   * f_j.x +
			   (rijzbar * rijybar)   * f_j.y +
			   (1+ rijzbar * rijzbar)* f_j.z  ) / mag_ijbar;

	    // deltaGim:
	      v_hydr.x += -2*r_i.z*r_j.z*(1/mag_ijbar3 - 3*(r_i.x - r_j.x)*(r_i.x - r_j.x)/mag_ijbar5)* f_j.x + 
			  6*(r_i.z*r_j.z*(r_i.x - r_j.x)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.y + 
			  4*(r_i.x - r_j.x)*(r_j.z/mag_ijbar3 - 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
	      
	      v_hydr.y += 6*(r_i.z*r_j.z*(r_i.x - r_j.x)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.x -
			  2*r_i.z*r_j.z*(1/mag_ijbar3 - 3*(r_i.y - r_j.y)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.y + 
			  4*(r_i.y - r_j.y)*(r_j.z/mag_ijbar3 - 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
	      
		
	      v_hydr.z += 2*(r_i.x - r_j.x)*(r_j.z/mag_ijbar3 + 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.x +
			  2*(r_i.y - r_j.y)*(r_j.z/mag_ijbar3 + 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.y +
			  4*r_i.z*r_j.z*(1/mag_ijbar3- 3*(r_i.z + r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
	    }
	  
	  
	}


      v_h[i].x = ( v_hydr.x ) *  (1./(8.*PI*eta));
      v_h[i].y = ( v_hydr.y ) *  (1./(8.*PI*eta));
      v_h[i].z = ( v_hydr.z ) *  (1./(8.*PI*eta));
//      v_h[i].x = (0.75 * v_hydr.x ) *  (hydr_radius/ gamma) * bond;
//      v_h[i].y = (0.75 * v_hydr.y ) *  (hydr_radius/ gamma) * bond;
//      v_h[i].z = (0.75 * v_hydr.z ) *  (hydr_radius/ gamma) * bond;
	
    }




}



/*  ************************************************************************ */

void     blake_2filaments( vec_s *r,
    vec_s *fa,
    vec_s *v_h,
    param_s parameters)
/*  ************************************************************************ */
{
  static float hydr_radius, eta, gamma, bond, massa, massa1, massa2;
  static char  init_flag = TRUE;
  static int n_mon1, n_mon2, n_fil;
  
  float rijx, rijy,rijz,rijxbar, rijybar,rijzbar;
  float mag_ij, mag_ijbar, mag_ijbar2, mag_ijbar3, mag_ijbar5;
  int  i,j,k; 
  vec_s r_i, r_j, f_i, f_j, v_hydr;

  if( init_flag )
    {
      eta            = parameters.viscosity;
      gamma          = parameters.gamma1*parameters.mass1;
      massa1          = parameters.mass1;
      massa2          = parameters.mass2;
      n_fil          = parameters.n_fil;
      n_mon1          = parameters.n_mon1;        
      n_mon2          = parameters.n_mon2;        
      init_flag  = FALSE;
    }



  /*    //printf("hydrodynamic radius is %f", hydr_radius); */
  
  for( i = 0; i < (n_mon1 + n_mon2); i++ )
    {

      v_hydr.x = 0;
      v_hydr.y = 0;
      v_hydr.z = 0;
      
      
      r_i.x      = r[i].x;
      r_i.y      = r[i].y;
      r_i.z      = r[i].z; 




      for( j = 0; j < (n_mon1 + n_mon2); j++)
  {
    
      if(j<n_mon1){
        massa = massa1;
      }
      else{
        massa = massa2;        
      }

    f_j.x      = massa * fa[j].x;
    f_j.y      = massa * fa[j].y;
    f_j.z      = massa * fa[j].z;
    
    r_j.x      = r[j].x;
    r_j.y      = r[j].y;
    r_j.z      = r[j].z;
  
    rijx     = r[i].x - r[j].x;
    rijy     = r[i].y - r[j].y;
    rijz     = r[i].z - r[j].z;
    rijxbar     = r[i].x - r[j].x;
    rijybar     = r[i].y - r[j].y;
    rijzbar     = r[i].z + r[j].z;

    mag_ij  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
    mag_ijbar2  = rijx*rijx + rijy*rijy + rijzbar*rijzbar;
    mag_ijbar  = sqrt(mag_ijbar2);
    mag_ijbar3 = mag_ijbar*mag_ijbar2;
    mag_ijbar5 = mag_ijbar3*mag_ijbar2;
    
    rijx    /= mag_ij;
    rijy    /= mag_ij;
    rijz    /= mag_ij;
    rijxbar  /= mag_ijbar;
    rijybar  /= mag_ijbar;
    rijzbar  /= mag_ijbar;
    
//    mag_ij  /= bond;
//    mag_ijbar  /= bond;
    
    
    
    if(j != i)
      {
      // OSEEN:
        v_hydr.x +=  ((1+ rijx * rijx)* f_j.x +
         (rijx * rijy)   * f_j.y +
         (rijx * rijz)   * f_j.z  ) / mag_ij;
        
        v_hydr.y +=  ((rijx* rijy)    * f_j.x +
         (1+ rijy * rijy)* f_j.y +
         (rijy * rijz)   * f_j.z  ) / mag_ij;
    
        v_hydr.z +=  ((rijz * rijx)   * f_j.x +
         (rijz * rijy)   * f_j.y +
         (1+ rijz * rijz)* f_j.z  ) / mag_ij;
         
      // OSEEN_BAR:
        v_hydr.x -= ((1+ rijxbar * rijxbar)* f_j.x +
         (rijxbar * rijybar)   * f_j.y +
         (rijxbar * rijzbar)   * f_j.z  ) / mag_ijbar;
        
        v_hydr.y -= ((rijxbar* rijybar)    * f_j.x +
         (1+ rijybar * rijybar)* f_j.y +
         (rijybar * rijzbar)   * f_j.z  ) / mag_ijbar;
    
        v_hydr.z -= ((rijzbar * rijxbar)   * f_j.x +
         (rijzbar * rijybar)   * f_j.y +
         (1+ rijzbar * rijzbar)* f_j.z  ) / mag_ijbar;

      // deltaGim:
        v_hydr.x += -2*r_i.z*r_j.z*(1/mag_ijbar3 - 3*(r_i.x - r_j.x)*(r_i.x - r_j.x)/mag_ijbar5)* f_j.x + 
        6*(r_i.z*r_j.z*(r_i.x - r_j.x)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.y + 
        4*(r_i.x - r_j.x)*(r_j.z/mag_ijbar3 - 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
        
        v_hydr.y += 6*(r_i.z*r_j.z*(r_i.x - r_j.x)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.x -
        2*r_i.z*r_j.z*(1/mag_ijbar3 - 3*(r_i.y - r_j.y)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.y + 
        4*(r_i.y - r_j.y)*(r_j.z/mag_ijbar3 - 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
        
    
        v_hydr.z += 2*(r_i.x - r_j.x)*(r_j.z/mag_ijbar3 + 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.x +
        2*(r_i.y - r_j.y)*(r_j.z/mag_ijbar3 + 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.y +
        4*r_i.z*r_j.z*(1/mag_ijbar3- 3*(r_i.z + r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
      }
    
    
  }

      v_h[i].x = ( v_hydr.x ) *  (1./(8.*PI*eta));
      v_h[i].y = ( v_hydr.y ) *  (1./(8.*PI*eta));
      v_h[i].z = ( v_hydr.z ) *  (1./(8.*PI*eta));
    
   // printf("part %i, vhyd = %f\n", i,v_h[i].z); 
    }




}


/* ########################################################################## */
float      boundary_conditions( vec_s *r,
                                vec_s *v,
                                vec_s *f,
                                param_s params,
                                int n_step,
                                int sstp )
/* ########################################################################## */
{
    static  vec_s *r_old;
    static  int   init_flag = TRUE;
    static  float b_l, w_length, freq, amp, l;
    static  int   n_mon, n_left, n_right, n_fil; 
    static  float gamma, hydr_radius;
    
    int   j,i;
    float v_new, t;
    float rijx, rijy, r0x, r0y, r0z, arg;
    float fac, x, dprod;
    float fx_tot, fy_tot, fz_tot;

/*ILl*/float aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom, denom1, num;
    vec_s v_per, v_par;

    if( init_flag )
    {
        n_mon    = params.n_mon;
	n_fil    = params.n_fil;
        b_l      = params.b_l;
        n_left   = params.n_left;
        n_right  = params.n_right;
        w_length = params.w_length;
        freq     = params.freq;
         amp      = params.amp;
        gamma    = params.gamma;
	hydr_radius = params.hydr_radius;
        l        = 1.0 - n_left*b_l - n_right*b_l;
        dx       = (float *)malloc( params.n_mon*sizeof( float ));
        dy       = (float *)malloc( params.n_mon*sizeof( float ));
        r_old    = (vec_s *)malloc( params.n_mon*sizeof( vec_s ));
        f_old    = (vec_s *)malloc( params.n_mon*sizeof( vec_s ));
        for( j = 0; j < params.n_mon; j++ )
	        {
           dx[j] = 0.0;
           dy[j] = 0.0;
        }
#ifdef EQUIL 
	freq      = 0.0;
#endif

        init_flag = FALSE;
    }


    t      = n_step*params.dt;
    w       = 0.0;


    for( i = n_left + 1; i < n_mon - n_right - 1; i++ )
    {
       x         = (i - n_left)*b_l/l;
#ifdef SQUARE 
       fac     = cos( w_length*2.0*PI*x - freq*t);
       if( fac > 0.0 )
          fac =  1.0;
       else
          fac = -1.0;
        
       fac  *= amp;
#endif

#ifdef STARTMODE 

      /*   QUI! metti if cycle <relax fac = amp cos(wl 2 Pi x)*/
       //nb nonfzia se amp aumenta con Sp!!!
       if(sstp < 20){
       fac       =  amp*cos( w_length*2.0*PI*x);
       }
       else{
       fac       =  amp*cos( w_length*2.0*PI*x - freq*t);
       }

#else

       fac       =  amp*cos( w_length*2.0*PI*x - freq*t);


#endif
       //       t_ext[i]  = fac;

 

	   r0x       =  (r[i].x + r[i+1].x)/2.0;
           r0y       =  (r[i].y + r[i+1].y)/2.0;
	   r0z       =  (r[i].z + r[i+1].z)/2.0;
	   
	   aa1 =  (r[i+1].x - r0x);
	   bb1 =  (r[i+1].y - r0y);
	   cc1 =  (r[i+1].z - r0z);
	   
	   aa2 =  (r[i+1].y - r0y)*(r0z - r[i].z) - (r[i+1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i+1].z - r0z)*(r0x - r[i].x) - (r[i+1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i+1].x - r0x)*(r0y - r[i].y) - (r[i+1].y - r0y)*(r0x - r[i].x);

	   denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }
	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
	   num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);

	   
       f[i+1].x += num*A1*fac/denom;
       f[i+1].y += num*fac/denom;
	   f[i+1].z += num*A2*fac/denom;
	   
	   //printf("11\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);
	   
	   aa1 =  (r0x - r[i].x);
	   bb1 =  (r0y - r[i].y);
	   cc1 =  (r0z - r[i].z);
	   
	   aa2 =  (r[i+1].y - r0y)*(r0z - r[i].z) - (r[i+1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i+1].z - r0z)*(r0x - r[i].x) - (r[i+1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i+1].x - r0x)*(r0y - r[i].y) - (r[i+1].y - r0y)*(r0x - r[i].x);


       denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }

	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
	   num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);

	   f[i].x += num*A1*fac/denom;
       f[i].y += num*fac/denom;
	   f[i].z += num*A2*fac/denom;

	  // //printf("22\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);



       r0x       =  (r[i].x + r[i-1].x)/2.0;
       r0y       =  (r[i].y + r[i-1].y)/2.0;
	   r0z       =  (r[i].z + r[i-1].z)/2.0;
       

	   aa1 =  -(r[i-1].x - r0x);
	   bb1 =  -(r[i-1].y - r0y);
	   cc1 =  -(r[i-1].z - r0z);
	   
	   aa2 =  (r[i-1].y - r0y)*(r0z - r[i].z) - (r[i-1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i-1].z - r0z)*(r0x - r[i].x) - (r[i-1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i-1].x - r0x)*(r0y - r[i].y) - (r[i-1].y - r0y)*(r0x - r[i].x);

	   denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }
	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
	 num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);
  
	   
       f[i-1].x += num*A1*fac/denom;
       f[i-1].y += num*fac/denom;
	   f[i-1].z += num*A2*fac/denom;

	  // //printf("33\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);


	   aa1 =  (r[i].x - r0x);
	   bb1 =  (r[i].y - r0y);
	   cc1 =  (r[i].z - r0z);
	   
	   aa2 =  (r[i-1].y - r0y)*(r0z - r[i].z) - (r[i-1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i-1].z - r0z)*(r0x - r[i].x) - (r[i-1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i-1].x - r0x)*(r0y - r[i].y) - (r[i-1].y - r0y)*(r0x - r[i].x);

	    denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }
	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);
       
	   f[i].x += num*A1*fac/denom;
       f[i].y += num*fac/denom;
	   f[i].z += num*A2*fac/denom;
	   
	   //printf("44\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);

    
	



/*
       r0x       =  (r[i].x + r[i+1].x)/2.0;
       r0y       =  (r[i].y + r[i+1].y)/2.0;
       rijx      = +(r[i+1].x - r0x);
       rijy      = -(r[i+1].y - r0y);
       f[i+1].x += rijy*fac;
       f[i+1].y += rijx*fac;

       rijx      = +(r[i].x - r0x);
       rijy      = -(r[i].y - r0y);
       f[i].x += rijy*fac;
       f[i].y += rijx*fac;

       r0x       =  (r[i].x + r[i-1].x)/2.0;
       r0y       =  (r[i].y + r[i-1].y)/2.0;
       rijx      =  -(r[i-1].x - r0x);
       rijy      =  +(r[i-1].y - r0y);

       f[i-1].x += rijy*fac;
       f[i-1].y += rijx*fac;
       rijx      =  -(r[i].x - r0x);
       rijy      =  +(r[i].y - r0y);
       f[i].x   += rijy*fac;
       f[i].y   += rijx*fac;*/
	   
	   
    }
    if( n_step > 0 )
    {
      fx_tot = 0.0;
      fy_tot = 0.0;
      
      for( i = 0; i < n_mon; i++ )
	{
	  fx_tot += f[i].x;
	  fy_tot += f[i].y;
	  
       }
      fx_tot /= n_mon;
      fy_tot /= n_mon;
      
 

      for( i = 0; i < n_mon; i++ )
       {
          dx[i]  = r[i].x - r_old[i].x;
          dy[i]  = r[i].y - r_old[i].y;
          w  += params.mass*(dx[i]*f_old[i].x + dy[i]*f_old[i].y);
          r_old[i].x = r[i].x;
          r_old[i].y = r[i].y;
          f[i].x -= fx_tot;
	      f[i].y -= fy_tot;
          f_old[i].x = f[i].x;
          f_old[i].y = f[i].y; 

       }
    }
 

  for( i = n_left + n_mon + 1; i < n_fil*n_mon - n_right - 1; i++ )
    {
       x         = (i - n_left - n_mon)*b_l/l;
#ifdef SQUARE 
       fac     = cos( w_length*2.0*PI*x - freq*t);
       if( fac > 0.0 )
          fac =  1.0;
       else
          fac = -1.0;
        
       fac  *= amp;
#endif

#ifdef STARTMODE 

      /*   QUI! metti if cycle <relax fac = amp cos(wl 2 Pi x)*/
       //nb nonfzia se amp aumenta con Sp!!!
       if(sstp < 20){
       fac       =  amp*cos( w_length*2.0*PI*x);
       }
       else{
       fac       =  amp*cos( w_length*2.0*PI*x - freq*t);
       }

#else

       fac       =  amp*cos( w_length*2.0*PI*x - freq*t);


#endif
       //       t_ext[i]  = fac;

   /*    r0x       =  (r[i].x + r[i+1].x)/2.0;
       r0y       =  (r[i].y + r[i+1].y)/2.0;
       rijx      = +(r[i+1].x - r0x);
       rijy      = -(r[i+1].y - r0y);
       f[i+1].x += rijy*fac;
       f[i+1].y += rijx*fac;

       rijx      = +(r[i].x - r0x);
       rijy      = -(r[i].y - r0y);
       f[i].x += rijy*fac;
       f[i].y += rijx*fac;

       r0x       =  (r[i].x + r[i-1].x)/2.0;
       r0y       =  (r[i].y + r[i-1].y)/2.0;
       rijx      =  -(r[i-1].x - r0x);
       rijy      =  +(r[i-1].y - r0y);

       f[i-1].x += rijy*fac;
       f[i-1].y += rijx*fac;
       rijx      =  -(r[i].x - r0x);
       rijy      =  +(r[i].y - r0y);
       f[i].x   += rijy*fac;
       f[i].y   += rijx*fac;*/
	   
	  
	   r0x       =  (r[i].x + r[i+1].x)/2.0;
       r0y       =  (r[i].y + r[i+1].y)/2.0;
	   r0z       =  (r[i].z + r[i+1].z)/2.0;
	   
	   aa1 =  (r[i+1].x - r0x);
	   bb1 =  (r[i+1].y - r0y);
	   cc1 =  (r[i+1].z - r0z);
	   
	   aa2 =  (r[i+1].y - r0y)*(r0z - r[i].z) - (r[i+1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i+1].z - r0z)*(r0x - r[i].x) - (r[i+1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i+1].x - r0x)*(r0y - r[i].y) - (r[i+1].y - r0y)*(r0x - r[i].x);

	   denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }
	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
	   num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);

	   
       f[i+1].x += num*A1*fac/denom;
       f[i+1].y += num*fac/denom;
	   f[i+1].z += num*A2*fac/denom;
	   
	  // //printf("2\t11\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);
	   
	   aa1 =  (r0x - r[i].x);
	   bb1 =  (r0y - r[i].y);
	   cc1 =  (r0z - r[i].z);
	   
	   aa2 =  (r[i+1].y - r0y)*(r0z - r[i].z) - (r[i+1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i+1].z - r0z)*(r0x - r[i].x) - (r[i+1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i+1].x - r0x)*(r0y - r[i].y) - (r[i+1].y - r0y)*(r0x - r[i].x);


       denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }

	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
	   num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);

	   f[i].x += num*A1*fac/denom;
       f[i].y += num*fac/denom;
	   f[i].z += num*A2*fac/denom;

	   //printf("2\t22\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);



       r0x       =  (r[i].x + r[i-1].x)/2.0;
       r0y       =  (r[i].y + r[i-1].y)/2.0;
	   r0z       =  (r[i].z + r[i-1].z)/2.0;
       

	   aa1 =  -(r[i-1].x - r0x);
	   bb1 =  -(r[i-1].y - r0y);
	   cc1 =  -(r[i-1].z - r0z);
	   
	   aa2 =  (r[i-1].y - r0y)*(r0z - r[i].z) - (r[i-1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i-1].z - r0z)*(r0x - r[i].x) - (r[i-1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i-1].x - r0x)*(r0y - r[i].y) - (r[i-1].y - r0y)*(r0x - r[i].x);

	   denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }
	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
	   num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);

	   
       f[i-1].x += num*A1*fac/denom;
       f[i-1].y += num*fac/denom;
	   f[i-1].z += num*A2*fac/denom;

	  // //printf("2\t33\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);


	   aa1 =  (r[i].x - r0x);
	   bb1 =  (r[i].y - r0y);
	   cc1 =  (r[i].z - r0z);
	   
	   aa2 =  (r[i-1].y - r0y)*(r0z - r[i].z) - (r[i-1].z - r0z)*(r0y - r[i].y);
	   bb2 =  (r[i-1].z - r0z)*(r0x - r[i].x) - (r[i-1].x - r0x)*(r0z - r[i].z);
	   cc2 =  (r[i-1].x - r0x)*(r0y - r[i].y) - (r[i-1].y - r0y)*(r0x - r[i].x);

	    denom1 = aa2*cc1-aa1*cc2;

	   if(denom1==0.0)
		{A1=0.0;
		A2=0.0;}
	   else
		{

		A1 = (bb1*cc2-bb2*cc1)/(aa2*cc1-aa1*cc2);
		A2 = (bb2*aa1-bb1*aa2)/(aa2*cc1-aa1*cc2);
	   }
	   
	   denom = sqrt(1.0+A1*A1+A2*A2);
num = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);

       
	   f[i].x += num*A1*fac/denom;
       f[i].y += num*fac/denom;
	   f[i].z += num*A2*fac/denom;
	   
	   //printf("2\t44\taa1=%e aa2=%e bb1=%e bb2=%e cc1=%e cc2=%e A1=%e A2=%e denom=%e\n", aa1, aa2, bb1, bb2, cc1, cc2, A1, A2, denom);

    
 
	   
    }
  if( n_step > 0 )
    {
      fx_tot = 0.0;
      fy_tot = 0.0;
      
      for( i = n_mon; i < n_fil*n_mon; i++ )
	{
	  fx_tot += f[i].x;
	  fy_tot += f[i].y;
	  
       }
      fx_tot /= n_mon;
      fy_tot /= n_mon;

    for( i = n_mon; i < n_fil*n_mon; i++ )
       {
          dx[i]  = r[i].x - r_old[i].x;
          dy[i]  = r[i].y - r_old[i].y;
          w  += params.mass*(dx[i]*f_old[i].x + dy[i]*f_old[i].y);
          r_old[i].x = r[i].x;
          r_old[i].y = r[i].y;
	  f[i].x -= fx_tot;
	  f[i].y -= fy_tot;
          f_old[i].x = f[i].x;
          f_old[i].y = f[i].y;
	  
       }


    }
    return( w );

}     
 

/* ########################################################################## */
void     output_data(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *f,
		       param_s params,
                       int     step,
                       float   *fper_ave,
                       float   *fpar_ave       )
/* ########################################################################## */
{
    FILE  *out_file_1, *out_file_2, *out_file_3;
    char  file_name[MAX_BUFF];
    int   i,j, n_mon;
    float rijx, rijy, rijz, rij, t, frac, dprod;
    float f_parx_tot, f_pary_tot, f_parz_tot;
    float f_perx_tot, f_pery_tot, f_perz_tot;
    float hydr_radius, gamma, mag;

    vec_s f_par, f_per;

    
    t    = step*params.dt;
    hydr_radius = params.hydr_radius;
    gamma = params.gamma;
    frac = t*params.freq/(2.0*PI); 
    frac = step/params.n_cycle;
    n_mon = params.n_mon;

    sprintf( file_name,"%s%s%.2f", params.dir_name, "out_", frac );
    out_file_1 = fopen( file_name, "w" );


  for(j=0; j<params.n_fil;j++){   
    //j=0;      
    for( i = 0; i < n_mon; i++ )
	{
	fprintf( out_file_1, "%.10e %.10e %.10e %.10e %.10e %.10e \n", 
		 r[i+j*n_mon].x, r[i+j*n_mon].y, r[i+j*n_mon].z, v[i+j*n_mon].x, v[i+j*n_mon].y, v[i+j*n_mon].z );
	}
      //fprintf( out_file_1, "\n");
    }
    fclose( out_file_1 );

} 


/* ########################################################################## */
void     output_data_2filaments(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *f,
           param_s params,
                       int     step,
                       float   *fper_ave,
                       float   *fpar_ave       )
/* ########################################################################## */
{
    FILE  *out_file_1, *out_file_2, *out_file_3;
    char  file_name[MAX_BUFF];
    int   i,j, n_mon1, n_mon2;
    float rijx, rijy, rijz, rij, t, frac, dprod;
    float f_parx_tot, f_pary_tot, f_parz_tot;
    float f_perx_tot, f_pery_tot, f_perz_tot;
    float hydr_radius, gamma, mag;

    vec_s f_par, f_per;

    
    t    = step*params.dt;
    frac = t*params.freq/(2.0*PI); 
    frac = step/params.n_cycle;
    n_mon1 = params.n_mon1;
    n_mon2 = params.n_mon2;

    sprintf( file_name,"%s%s%.2f", params.dir_name, "out_", frac );
    out_file_1 = fopen( file_name, "w" );

    for( i = 0; i < (n_mon1+n_mon2); i++ )
  {
    fprintf( out_file_1, "%.10e %.10e %.10e %.10e %.10e %.10e \n", 
       r[i].x, r[i].y, r[i].z, v[i].x, v[i].y, v[i].z );
  }
      //fprintf( out_file_1, "\n");
    
    fclose( out_file_1 );

} 


void print_force(vec_s *force, param_s par, int n_mon, char *filename)
{
   int i;
   static int num=0;
   FILE *f_ptr;
   char fname[100];
   char command[100];
   
  /*      sprintf(command,"mkdir -p %s%s",par.dir_name, "Force/"); */
/*        system(command); */
      
      sprintf(fname,"%s%s%s_%d",par.dir_name, "Force/", filename, num);

      f_ptr = fopen(fname, "w");
      for(i=0;i<=n_mon-1;i++){
         fprintf(f_ptr,"%d %.16g %.16g \n", i, force[i].x, force[i].y);
      }
      fclose(f_ptr);
   num++;
 
}

