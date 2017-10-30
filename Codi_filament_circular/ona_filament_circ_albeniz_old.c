#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <ran.h>*/
#include "useful.h"
#include "../../../0_Codi/Funcions/os3many_filament.h"
#include "../../../0_Codi/Funcions/vectors.h"
#include "../../../0_Codi/Funcions/vectors.c"
#include "../../../0_Codi/Funcions/funcions_filament.c"
#define BFIELD
//#define BUCKLING_INDUIT
//#define BUCKLING
//#define TRANS
//#define CONSTANT_FORCE
//#define HEAD
int main( int argc, char *argv[] )
{
    vec_s *r, *v;

    param_s parameters;
    configure_sys( argc, argv, &parameters, &r, &v);
    simulate( parameters, r, v );

    return(1);
}

/* ######################################################################### */
void  configure_sys( int      argc,
                     char    *argv[],
                     param_s *parameters,
                     vec_s  **r,
                     vec_s  **v  )
/* ######################################################################### */
{
    FILE  *param_file, *info_file, *config_file;

    char   buffer[MAX_BUFF], file_name[MAX_BUFF], dir_name[MAX_BUFF];
    float  sp_numb4, tinercia, Fmax;
    int    n_mon,n_fil,i;
    float  rijx, rijy, rijz, B;

    seed = -789219;

    if( argc == 1 )
    {
             // printf("\nNo parameter file name\n" );
       exit( 1 );
    }
    param_file = fopen( argv[1], "r" );
    if( param_file == NULL )
    {
        // printf("\nParameter file not found\n");
         exit( 1 );
    }
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%s", file_name );
    config_file = fopen( file_name, "r" );
    if( config_file == NULL )
    {
        // printf("\nConfiguration file %s not found\n", file_name );
         fflush( stdout );
         exit( 1 );
    }

    sprintf( parameters->outfile_name,"%s%s", file_name, "_new" );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%s", dir_name );
    sprintf( file_name,"%s%s", dir_name, "info");
    sprintf( parameters->dir_name,"%s", dir_name );
    info_file = fopen( file_name, "w" );
    if( info_file == NULL )
    {
        // printf("\nCannot access directory for output %s\n", file_name );
         fflush( stdout );
         exit( 1 );
    }
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->n_steps) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->dt) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->k) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->mass) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->viscosity) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->gamma) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->hydr_radius) );
    fgets( buffer, MAX_BUFF, param_file );

    sscanf( buffer, "%le", &(parameters->tolerance) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->n_cycle) );
    fgets( buffer, MAX_BUFF, param_file );

    sscanf( buffer, "%le", &(parameters->Fx) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->Fz) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->Bx) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->Bz) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->freq) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->chi) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->m2) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->mu0) );
    fgets( buffer, MAX_BUFF, param_file );

    sscanf( buffer, "%d", &(parameters->oseen) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->blake) );

    
    /* FM: Sperm Number computation. N.B. L=1 ! */
    sp_numb4 = 2*PI*(parameters->gamma)/(parameters->k);
    sp_numb4 /= ((parameters->n_cycle)*(parameters->dt));
    parameters->Sperm_Number = pow(sp_numb4, 0.25);
    
    

    n_mon = 0;
    n_fil = 1;
    parameters->n_fil = n_fil;
  
    while( fgets( buffer, MAX_BUFF, config_file ) != NULL )
    {
         n_mon++;
    }
   
    rewind(config_file);

    fflush( stdout );
    
    /*   posiz, velocita' di tutti i bead */  

    (*r) = (vec_s *)malloc(  n_mon*n_fil*sizeof( vec_s ));
    (*v) = (vec_s *)malloc(  n_mon*n_fil*sizeof( vec_s ));
       
    for(i=0;i<n_mon;i++){     
      fgets(  buffer, MAX_BUFF, config_file );
      sscanf( buffer,"%le%le%le%le%le%le", &((*r)[i].x),
                                              &((*r)[i].y),
                                              &((*r)[i].z),
                                              &((*v)[i].x),
                                              &((*v)[i].y),
                                              &((*v)[i].z)    );    
	

    }
    fclose( config_file );

    
    //// printf("n_mon %d, n_fil %d\n",n_mon, n_fil);

    rijx = (*r)[0].x - (*r)[1].x;
    rijy = (*r)[0].y - (*r)[1].y;
    rijz = (*r)[0].z - (*r)[1].z;
    parameters->b_l         = sqrt(rijx*rijx + rijy*rijy + rijz*rijz );
    parameters->n_mon       = n_mon;
    parameters->mass       /= (float)n_mon;

    B = parameters->Fz*parameters->b_l/parameters->k;
    parameters->gamma      /= parameters->mass*n_mon;
    parameters->k          /= parameters->b_l*parameters->b_l;
    parameters->k          /= (float)(n_mon - 2) ;
//    parameters->freq        = 2.0*PI/(float)parameters->n_cycle;
//    parameters->freq       /= parameters->dt;
//    parameters->freq =   0.1*parameters->gamma;
 
    fclose(info_file);
}

/* ########################################################################## */

void  simulate(     param_s  params,
                     vec_s  *r,
                     vec_s  *v)
/* ########################################################################## */
{
    vec_s *fc_1, *fc_2, *fc_3;
    vec_s *f_bend, *fc_1a, *fc_2a,*fc_3a, *f_dip;

    int    i, j, k, l,m, n_mon, step, n_ang, n_fil, nmig;
    float  *mags, *mags1;
    float  u, ke, kk;
    float  rijx, rijy, rijz;
    float  fpar_ave, fper_ave;
    float  cycle_time;
    float  sc1;
    vec_s  r0f, m2, vec1, vec2, Bfield;
    vec_s *rij, *rij1, *f,*fa,*v1, *r1, *rm, *rcm;

    FILE *conf = fopen("conf", "w");
   // FILE *flog_fa;
    /* Fab */ float freq, amplitude;

    /* Fab */ float velocity, velocityx,velocityy, velocityz, velx, vely,velz;
    /* Fab */ int   count;
    /* Fab */ float Ext_force;
    /* Fab */ float Sperm_Number;
    /* Fab */ float x_old, y_old, z_old;
    /* Fab */ char filename[100];
    /* Fab */ FILE *file_ptr, *logfile;

    /* MCL */ float  loc_dr, fatt, fatt2, b_l;
    /* MCL */ float  dr0x, dr0y,drijx,drijy,drijz;
    /* MCL */ float dist, Torque;

    FILE *frcm = fopen("RCM.dat", "w");

    /* Fab */ vec_s constant_force;
    /* Fab */ constant_force.x = 0.0;
    /* Fab */ constant_force.y = 0.0;
    /* Fab */ constant_force.z = 0.0;


       
    /* Fab */ freq = params.freq;

    
    dist  = params.dist;
    n_fil = params.n_fil;
    n_mon = params.n_mon;   
    dist = 0.5;
    
/* Forze  a = di tutti  */    

    f_bend   = (vec_s *)malloc( n_mon*sizeof( vec_s ));
    f_dip   = (vec_s *)malloc(sizeof( vec_s ));

    fc_1a    = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));  
    fc_1     = (vec_s *)malloc( n_mon*sizeof( vec_s ));

    fc_2a    = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));
    fc_2     = (vec_s *)malloc( n_mon*sizeof( vec_s ));
 
    fc_3a    = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));
    fc_3     = (vec_s *)malloc( n_mon*sizeof( vec_s ));    
    

    rij   = (vec_s *)malloc((n_mon + 1)*n_fil*sizeof( vec_s ));
    rij1  = (vec_s *)malloc((n_mon + 1)*sizeof( vec_s ));

    fa    = (vec_s *)malloc(n_mon*n_fil*sizeof( vec_s ));
    f     = (vec_s *)malloc(n_mon*sizeof( vec_s ));

    v1    = (vec_s *)malloc(n_mon*sizeof( vec_s ));
    r1    = (vec_s *)malloc(n_mon*sizeof( vec_s ));

    mags  = (float *)malloc((n_mon+1)*n_fil*sizeof( float ));
    mags1 = (float *)malloc((n_mon+1)*sizeof( float ));
  
    rm    = (vec_s *)malloc(n_mon*n_fil*sizeof( vec_s ));
    rcm = (vec_s *)malloc(n_fil*sizeof( vec_s ));


  //  sprintf(filename,"logfile.log",params.dir_name);
    logfile = fopen("logfile.log","w");


    /*  generates config  n_fil filaments */
    for(i=0; i<n_fil; i++){
      for(j=0; j<n_mon; j++){
	r[i*n_mon+j].x = r[j].x;
	r[i*n_mon+j].y = r[j].y;
	r[i*n_mon+j].z = r[j].z + i* dist;
	v[i*n_mon+j].x = v[j].x;
	v[i*n_mon+j].y = v[j].y;
	v[i*n_mon+j].z = v[j].z;
      }
    }

    
    
    cycle_time  = params.n_cycle*params.dt;

    for( i = 1; i < n_mon*n_fil; i++ )
    {
       rijx     = r[i-1].x - r[i].x;
       rijy     = r[i-1].y - r[i].y;
       rijz     = r[i-1].z - r[i].z;
       if(i%n_mon==0){
	 rijx     = 0;
	 rijy     = 0;
	 rijz     = 0;  
	 } 
     
       mags[i]  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
       rij[i].x = rijx;
       rij[i].y = rijy;
       rij[i].z = rijz;
    }

// AFEGIT PER FILAMENT CIRCULAR ----------
    rij[n_mon].x = r[n_mon-1].x - r[0].x;
    rij[n_mon].y = r[n_mon-1].y - r[0].y;
    rij[n_mon].z = r[n_mon-1].z - r[0].z;
    mags[n_mon] = sqrt(rij[n_mon].x*rij[n_mon].x + rij[n_mon].y*rij[n_mon].y + rij[n_mon].z*rij[n_mon].z);
   // printf("mags[n_mon] = %f\n",mags[n_mon]);
// ----------------------------------------



    for( i = 0; i < n_mon*n_fil; i++ )
    {
       fa[i].x = 0.0;
       fa[i].y = 0.0;
       fa[i].z = 0.0;

       fc_1a[i].x = 0.0;
       fc_1a[i].y = 0.0;
       fc_1a[i].z = 0.0;
       fc_2a[i].x = 0.0;
       fc_2a[i].y = 0.0;
       fc_2a[i].z = 0.0;
       fc_3a[i].x = 0.0;
       fc_3a[i].y = 0.0;
       fc_3a[i].z = 0.0;
    }

    for(m=0; m<n_fil; m++){
      for( i = 0; i < n_mon+1; i++ ){
    rij1[i].x = rij[i+m*n_mon].x;
    rij1[i].y = rij[i+m*n_mon].y;
    rij1[i].z = rij[i+m*n_mon].z;


    mags1[i] = mags[i+m*n_mon];

      }
      for( i = 0; i < n_mon; i++ ){

    f[i].x = fa[i+m*n_mon].x;
    f[i].y = fa[i+m*n_mon].y;
    f[i].z = fa[i+m*n_mon].z;
      }
    // printf("mags[n_mon] = %f\n",mags1[n_mon]);

      u    = bending_forces_circ( rij1, mags1, f, params );
    // printf("mags[n_mon] = %f\n",mags1[n_mon]);
      for( i = 0; i < n_mon; i++ ){
	fa[i+m*n_mon].x = f[i].x;
	fa[i+m*n_mon].y = f[i].y;
	fa[i+m*n_mon].z = f[i].z;
	
      }
    }


#ifdef BFIELD
    /* Fab */
    /* NO INTERACTION BETWEEN FILAMENTS !! */

        f_dip[0].x = 0.0;
        f_dip[0].y = 0.0;
        f_dip[0].z = 0.0;
        nmig = 15;
        Bfield.x = params.Bx*sin(freq*0.*params.dt);
        Bfield.z = params.Bz*cos(freq*0.*params.dt);
        Bfield.y = 0.0;
       // Bfield.z = 0.0;

        r0f.x = r[0].x - r[nmig].x;
        r0f.y = r[0].y - r[nmig].y;
        r0f.z = r[0].z - r[nmig].z;

/*        vec1 = vecAddition(rij[nmig], rij[nmig+1]);
        vec1 = scalarProduct(vec1, 0.5);
        sc1 = 1.0/sqrt(dotProduct(vec1,vec1));
        vec1 = scalarProduct(vec1, sc1);
        vec2.x = 0.;
        vec2.y = 0.;
        vec2.z = 1.0;
        m2 = crossProduct(vec1, vec2);
*/
        vec1.x = (r[nmig].x-r[nmig-1].x) + (r[nmig].x-r[nmig+1].x);
        vec1.y = (r[nmig].y-r[nmig-1].y) + (r[nmig].y-r[nmig+1].y);
        vec1.z = (r[nmig].z-r[nmig-1].z) + (r[nmig].z-r[nmig+1].z);
        sc1 = 1./sqrt(dotProduct(vec1,vec1));
        m2 = scalarProduct(vec1, sc1);
        m2 = scalarProduct(m2, params.m2);

        dipolar_forces(r0f, m2, f_dip, params, Bfield);
//        printf("%i\t%le\t%le\t%le\n", 0, m2.x, m2.y,m2.z);
        fa[0].x += f_dip[0].x/params.mass;
        fa[0].y += f_dip[0].y/params.mass;
        fa[0].z += f_dip[0].z/params.mass;
        fa[nmig].x -= f_dip[0].x/params.mass;
        fa[nmig].y -= f_dip[0].y/params.mass;
        fa[nmig].z -= f_dip[0].z/params.mass;

        Torque   =  (Bfield.x*m2.z - Bfield.z*m2.x)/params.mass; // Suposo que el camp magnetic es al pla XZ 
         drijx      =  -(r[nmig].x - r[0].x)/2.;
         drijz      =  +(r[nmig].z - r[0].z)/2.;
         
         fa[nmig].x += 0.5*drijz*Torque/(drijx*drijx + drijz*drijz);
         fa[nmig].z += 0.5*drijx*Torque/(drijx*drijx + drijz*drijz);
         drijx      =  -(r[0].x - r[nmig].x)/2.;
         drijz      =  +(r[0].z - r[nmig].z)/2.;
         fa[0].x   += 0.5*drijz*Torque/(drijx*drijx + drijz*drijz);
         fa[0].z   += 0.5*drijz*Torque/(drijx*drijx + drijz*drijz);
    
#endif



#ifdef BUCKLING_INDUIT
    /* Fab */
    /* Add constant Force on head of filament */

    fa[0].x += (params.Fx/params.mass)*sin(freq*0.0*params.dt);
    fa[n_mon-1].x -= (params.Fx/params.mass)*sin(freq*0.0*params.dt);
    nmig = (int)(n_mon/2.);
    if (sin(freq*0.0*params.dt) > 0 ){
      fa[nmig].z += (params.Fz/params.mass)*sin(freq*0.0*params.dt);
    }
    else{
      fa[nmig].z += 0.0;
    }
#endif
	     
	     /*  inizio ciclo temporale */

////////////////////////////////////////////////////////////////////////////////////////////
//COMENÇA EL CICLE TEMPORAL!!!!!!!!!!!!!!!
////////////////////////////////////////////////////////////////////////////////////////////

   for( j = 0; j < params.n_steps; j++ )
    {
       n_ang    = 0;
     
       /* Fab */
       velocity = 0.0;
       velocityx = 0.0;
       velocityy = 0.0;
       velocityz = 0.0;
       count    = 0;
    for( m = 0; m < n_fil; m++){
       for( i = 0; i < n_mon; i++ ){
	  rm[i + m*n_mon].x=0.0;
	  rm[i + m*n_mon].y=0.0;
	  rm[i + m*n_mon].z=0.0;}
      }

    for( i = 0; i < n_fil; i++ ){
	  rcm[i].x=0.0;
	  rcm[i].y=0.0;
	  rcm[i].z=0.0;}

    for( k = 0; k < params.n_cycle; k++ )
	 {
	   
	   /* FM: to compute the movement at each time step */
	   x_old = r[0].x;
	   y_old = r[0].y;
	   z_old = r[0].z;
    

          step = j*params.n_cycle + k;
	  
	  //fprintf(flog_fa,"  a: %i\t%le",step, fa[15].y);   

          fflush( stdout );

	  for(m=0; m<n_fil; m++){
	    for( i = 0; i < n_mon+1; i++ ){
	      rij1[i].x = rij[i+m*n_mon].x;
	      rij1[i].y = rij[i+m*n_mon].y;
	      rij1[i].z = rij[i+m*n_mon].z;
	 
	   
	      mags1[i] = mags[i+m*n_mon];
	
	    }

        for( i = 0; i < n_mon; i++ ){
        
          v1[i].x = v[i+m*n_mon].x;
          v1[i].y = v[i+m*n_mon].y;
          v1[i].z = v[i+m*n_mon].z;
    
        }
	     
    // printf("mags[n_mon] = %f\n",mags1[n_mon]);

	  // constrain_velocities( rij1, mags1, v1,fc_1, params );
	     constrain_velocities_circ( rij1, mags1, v1,fc_1, params );
	
	    for( i = 0; i < n_mon; i++ ){
	      v[i+m*n_mon].x = v1[i].x;
	      v[i+m*n_mon].y = v1[i].y;
	      v[i+m*n_mon].z = v1[i].z;
	     
	

	      fc_1a[i+m*n_mon].x   = fc_1[i].x;
	      fc_1a[i+m*n_mon].y   = fc_1[i].y;
	      fc_1a[i+m*n_mon].z   = fc_1[i].z;

	    }
	  }
	   
	
  
	  //fprintf(flog_fa,"  b: %le",fa[15].y);
          ke = verlet_pt1( r, v, fa, params );


	
	  for(m=0; m<n_fil; m++){
	    for( i = 0; i < n_mon+1; i++ ){
	      rij1[i].x = rij[i+m*n_mon].x;
	      rij1[i].y = rij[i+m*n_mon].y;
	      rij1[i].z = rij[i+m*n_mon].z;
	      
	      mags1[i] = mags[i+m*n_mon];

	
	    }    
        for( i = 0; i < n_mon; i++ ){
          v1[i].x = v[i+m*n_mon].x;
          v1[i].y = v[i+m*n_mon].y;
          v1[i].z = v[i+m*n_mon].z;

    
          r1[i].x = r[i+m*n_mon].x;
          r1[i].y = r[i+m*n_mon].y;
          r1[i].z = r[i+m*n_mon].z;
      }
        constrain_positions_circ( r1, rij1, v1, fc_2, params );
        //constrain_positions( r1, rij1, v1, fc_2, params );
	    
	    for( i = 0; i < n_mon; i++ ){
	      
	    v[i+m*n_mon].x = v1[i].x;
	    v[i+m*n_mon].y = v1[i].y;
	    v[i+m*n_mon].z = v1[i].z;
	    
	   
	    r[i+m*n_mon].x = r1[i].x;
	    r[i+m*n_mon].y = r1[i].y;
	    r[i+m*n_mon].z = r1[i].z;


	

	    fc_2a[i+m*n_mon].x   = fc_2[i].x;
	    fc_2a[i+m*n_mon].y   = fc_2[i].y;
	    fc_2a[i+m*n_mon].z   = fc_2[i].z;
	    }
	  }

	   
	  /*	  if( step % params.out_freq == 0 )
	    {
              output_data( r, v, fa, params, step, &fper_ave, &fpar_ave );
              fflush( stdout );
	    }
	  */

	  for( i = 1; i < n_mon*n_fil; i++ )
	    {
	      rijx     = r[i-1].x - r[i].x;
	      rijy     = r[i-1].y - r[i].y;
	      rijz     = r[i-1].z - r[i].z;
	      if(i%n_mon==0){
		rijx     = 0;
		rijy     = 0;
		rijz     = 0;
	      }
	      mags[i]  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
	      rij[i].x = rijx;
	      rij[i].y = rijy;
	      rij[i].z = rijz;
	    }
    // AFEGIT PER FILAMENT CIRCULAR ----------
        rij[n_mon].x = r[n_mon-1].x - r[0].x;
        rij[n_mon].y = r[n_mon-1].y - r[0].y;
        rij[n_mon].z = r[n_mon-1].z - r[0].z;
        mags[n_mon] = sqrt(rij[n_mon].x*rij[n_mon].x + rij[n_mon].y*rij[n_mon].y + rij[n_mon].z*rij[n_mon].z);
    // ----------------------------------------


	

	

          for( i = 0; i < n_mon*n_fil; i++ )
          {
	
	    fa[i].x = 0.0;
	    fa[i].y = 0.0;
	    fa[i].z = 0.0;
	  }


#ifdef BFIELD
    /* Fab */
    /* NO INTERACTION BETWEEN FILAMENTS !! */

        f_dip[0].x = 0.0;
        f_dip[0].y = 0.0;
        f_dip[0].z = 0.0;
        nmig = 15;
        Bfield.x = params.Bx*sin(freq*step*params.dt);
        Bfield.y = 0.0;
        Bfield.z = params.Bz*cos(freq*step*params.dt);
        //Bfield.z = 0.0;

        r0f.x = r[0].x - r[nmig].x;
        r0f.y = r[0].y - r[nmig].y;
        r0f.z = r[0].z - r[nmig].z;

/*        vec1 = vecAddition(rij[nmig], rij[nmig+1]);
        vec1 = scalarProduct(vec1, 0.5);
        sc1 = 1.0/sqrt(dotProduct(vec1,vec1));
        vec1 = scalarProduct(vec1, sc1);
        vec2.x = 0.;
        vec2.y = 0.;
        vec2.z = 1.0;
        m2 = crossProduct(vec1, vec2);
*/


        vec1.x = (r[nmig].x-r[nmig-1].x) + (r[nmig].x-r[nmig+1].x);
        vec1.y = (r[nmig].y-r[nmig-1].y) + (r[nmig].y-r[nmig+1].y);
        vec1.z = (r[nmig].z-r[nmig-1].z) + (r[nmig].z-r[nmig+1].z);
        sc1 = 1./sqrt(dotProduct(vec1,vec1));
        m2 = scalarProduct(vec1, sc1);
        m2 = scalarProduct(m2, params.m2);

        dipolar_forces(r0f, m2, f_dip, params, Bfield);

        fa[0].x += f_dip[0].x/params.mass;
        fa[0].y += f_dip[0].y/params.mass;
        fa[0].z += f_dip[0].z/params.mass;
        fa[nmig].x -= f_dip[0].x/params.mass;
        fa[nmig].y -= f_dip[0].y/params.mass;
        fa[nmig].z -= f_dip[0].z/params.mass;

        Torque   =  (Bfield.x*m2.z - Bfield.z*m2.x)/params.mass; // Suposo que el camp magnetic es al pla XZ 
         drijx      =  -(r[nmig].x - r[0].x)/2.;
         drijz      =  +(r[nmig].z - r[0].z)/2.;
         
         fa[nmig].x += 0.5*drijz*Torque/(drijx*drijx + drijz*drijz);
         fa[nmig].z += 0.5*drijx*Torque/(drijx*drijx + drijz*drijz);
         drijx      =  -(r[0].x - r[nmig].x)/2.;
         drijz      =  +(r[0].z - r[nmig].z)/2.;
         fa[0].x   += 0.5*drijz*Torque/(drijx*drijx + drijz*drijz);
         fa[0].z   += 0.5*drijz*Torque/(drijx*drijx + drijz*drijz);
    

 
#endif



#ifdef BUCKLING_INDUIT
    /* Fab */
    /* Add constant Force on head of filament */
    fa[0].x += (params.Fx/params.mass)*sin(freq*step*params.dt);
    fa[n_mon-1].x -= (params.Fx/params.mass)*sin(freq*step*params.dt);
    nmig = (int)(n_mon/2.);
    nmig = 15;
    if (sin(freq*step*params.dt) > 0 ){
      fa[nmig].z += (params.Fz/params.mass)*sin(freq*step*params.dt);
      fa[0].z -= 0.5*(params.Fz/params.mass)*sin(freq*step*params.dt);
      fa[n_mon-1].z -= 0.5*(params.Fz/params.mass)*sin(freq*step*params.dt);
    }
    else{
      fa[nmig].z = 0.0;
    }
#endif
	  

	  
	  for(m=0; m<n_fil; m++){
	    for( i = 0; i < n_mon+1; i++ ){

	      	

	      rij1[i].x = rij[i+m*n_mon].x;
	      rij1[i].y = rij[i+m*n_mon].y;
	      rij1[i].z = rij[i+m*n_mon].z;


	      mags1[i] = mags[i+m*n_mon];

	
	    }
        for( i = 0; i < n_mon; i++ ){

          
          f_bend[i].x = 0.0;
          f_bend[i].y = 0.0;
          f_bend[i].z = 0.0;
    
        }
     
	    u  = bending_forces_circ( rij1, mags1, f_bend, params );
     

	    for( i = 0; i < n_mon; i++ ){
	      fa[i+m*n_mon].x += f_bend[i].x;
	      fa[i+m*n_mon].y += f_bend[i].y;
	      fa[i+m*n_mon].z += f_bend[i].z;	
      
	    }
	  }


      for(m=0; m<n_fil; m++){
        for( i = 0; i < n_mon+1; i++ ){
          rij1[i].x = rij[i+m*n_mon].x;
          rij1[i].y = rij[i+m*n_mon].y;
          rij1[i].z = rij[i+m*n_mon].z;
     
       
          mags1[i] = mags[i+m*n_mon];
    
        }

        for( i = 0; i < n_mon; i++ ){

          f[i].x = fa[i+m*n_mon].x;
          f[i].y = fa[i+m*n_mon].y;
          f[i].z = fa[i+m*n_mon].z;
    
        }
         

	   
         // constrain_forces(rij1, mags1, f, fc_3, params );
          constrain_forces_circ(rij1, mags1, f, fc_3, params );
	
	    for( i = 0; i < n_mon; i++ ){
	      fc_3a[i+m*n_mon].x   = fc_3[i].x;
	      fc_3a[i+m*n_mon].y   = fc_3[i].y;
	      fc_3a[i+m*n_mon].z   = fc_3[i].z;

	    }
	  }
	   
	
          for( i = 0; i < n_mon*n_fil; i++ )
          {
	    fa[i].x  += fc_3a[i].x;
	    fa[i].y  += fc_3a[i].y;
	    fa[i].z  += fc_3a[i].z;
	   
	      
          }       
	 // fprintf(flog_fa,"  f: %le\n",fa[15].y);
	  verlet_pt2( r, v, fa , params );

      /* Fab */

	  /* Mean velocity at every cycle*/

	  velx=0.0;
	  vely=0.0;
	  velz=0.0;

//printf("\n");
	  for( i = 0; i < n_mon; i++ ){
	  //printf("%e\t%e\t%e\n", v[i].x, v[i].y, v[i].z);
	      velx += v[i].x;
	      vely += v[i].y;	  
	      velz += v[i].z;}
//printf("\n");

	  velx /= n_mon;
	  vely /= n_mon;
	  velz /= n_mon;

	  velocityx += velx;
	  velocityy += vely;
	  velocityz += velz;


	  /*Position of the CM*/
	
  for(m=0; m<n_fil; m++){
	  for( i = 0; i < n_mon; i++ ){
	  rm[i+n_mon*m].x += r[i+n_mon*m].x;
	  rm[i+n_mon*m].y += r[i+n_mon*m].y;
	  rm[i+n_mon*m].z += r[i+n_mon*m].z;}
  }

	  
if(step%params.n_cycle==0){
        
    for( i = 0; i < n_fil; i++ ){
   
          rcm[i].x = 0.0;
          rcm[i].y = 0.0;
          rcm[i].z = 0.0;
    }  

    for( i = 0; i < n_fil; i++ ){
        for(m=0;m<n_mon; m++){
          rcm[i].x += r[m].x/n_mon;
          rcm[i].y+= r[m].y/n_mon;
          rcm[i].z += r[m].z/n_mon;}
    }

      fprintf(frcm,"%i\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", step, rcm[0].x, rcm[0].y,rcm[0].z, r[0].x, r[0].y,r[0].z, r[nmig].x, r[nmig].y,r[nmig].z, Bfield.x );

    for( i = 0; i < n_fil; i++ ){
   
          rcm[i].x = 0.0;
          rcm[i].y = 0.0;
          rcm[i].z = 0.0;
    }        

    }






    }
    
/*  fine ciclo temporale */

       
      
       
     /* Fab */ /* Accumulates average over cycles */
       velocityx /= params.n_cycle;
       velocityy /= params.n_cycle;
       velocityz /= params.n_cycle;
       velocity = sqrt(velocityx*velocityx+velocityy*velocityy+velocityz*velocityz);
	  
       	 //if(j==params.n_steps-1)   
       kk = sin(freq*step*params.dt);
  //  fprintf(logfile, "%i\t%le\t%le\t%le\t%le\n", j, Bfield.x, Bfield.y,Bfield.z, velocityx);  
	printf("%i\t%le\t%le\t%le\t%le\n", j, velocity, velocityx,velocityy, velocityz );
       
       velocityx = 0.0;
       velocityy = 0.0;
       velocityz = 0.0;
	

//       velocity = (r[n_mon-1].z - z_start)/cycle_time;
      for( m = 0; m < n_fil; m++){
       for( i = 0; i < n_mon; i++ ){
	 rm[i+n_mon*m].x /= params.n_cycle;
	 rm[i+n_mon*m].y /= params.n_cycle;
	 rm[i+n_mon*m].z /= params.n_cycle;}
       }
       
       output_data( rm, v, fa, params, step, &fper_ave, &fpar_ave );
       fflush( stdout );
	  





   }
       fclose(logfile);
 //  fclose(flog_fa);
/*  fine tutti cicli  */

fclose(frcm);



}



