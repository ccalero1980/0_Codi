#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <ran.h>*/
#include "useful.h"
#include "os3many.h"

#include "/media/sf_carloscalero2/Dropbox/ABIOMATER/Nedador_individual/Filament_magnetic/Codi/Funcions/funcions.c"


#define TRANS
#define CONSTANT_FORCE
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
              printf("\nNo parameter file name\n" );
       exit( 1 );
    }
    param_file = fopen( argv[1], "r" );
    if( param_file == NULL )
    {
         printf("\nParameter file not found\n");
         exit( 1 );
    }
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%s", file_name );
    config_file = fopen( file_name, "r" );
    if( config_file == NULL )
    {
         printf("\nConfiguration file %s not found\n", file_name );
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
         printf("\nCannot access directory for output %s\n", file_name );
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
    sscanf( buffer, "%le", &(parameters->gamma) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->hydr_radius) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->gamma_head) );
    fgets( buffer, MAX_BUFF, param_file );

    sscanf( buffer, "%le", &(parameters->phi_0) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->tolerance) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->amp) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->n_cycle) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->w_length) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->n_left) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->n_right) );

    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->out_freq) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->amplitude) );

    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->Fx) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->Fz) );
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

    
    // printf("n_mon %d, n_fil %d\n",n_mon, n_fil);

    rijx = (*r)[0].x - (*r)[1].x;
    rijy = (*r)[0].y - (*r)[1].y;
    rijz = (*r)[0].z - (*r)[1].z;
    parameters->b_l         = sqrt(rijx*rijx + rijy*rijy + rijz*rijz );
    parameters->n_mon       = n_mon;
    parameters->mass       /= (float)n_mon;
#ifdef HEAD
    parameters->amplitude /= parameters->mass;
#endif

    B = parameters->Fz*parameters->b_l/parameters->k;
    parameters->gamma      /= parameters->mass*n_mon;
    parameters->gamma_head /= parameters->mass;
    parameters->amp        *= parameters->k;
    parameters->amp        /= parameters->mass;
    parameters->amp        /= parameters->b_l*parameters->b_l;
    parameters->k          /= parameters->b_l*parameters->b_l;
    parameters->k          /= (float)(n_mon - 2) ;
    parameters->freq        = 2.0*PI/(float)parameters->n_cycle;
    parameters->freq       /= parameters->dt;


  


        
    fprintf( info_file, "Number of steps   = %d\n", parameters->n_steps  );
    fprintf( info_file, "Time-step         = %e\n", parameters->dt       );
    fprintf( info_file, "Bending modulus   = %e\n", parameters->k        );
    fprintf( info_file, "Equilibrium angle = %e\n", parameters->phi_0    ); 
    fprintf( info_file, "Bond length       = %e\n", parameters->b_l      );
    fprintf( info_file, "Tolerance         = %e\n", parameters->tolerance); 
    fprintf( info_file, "hydr radius       = %e\n", parameters->hydr_radius *parameters->b_l ); 
    fprintf( info_file, "Friction gamma    = %e\n", parameters->gamma    );
    fprintf( info_file, "Frequency         = %e\n", parameters->freq     );
    fprintf( info_file, "Peak torque       = %e\n", parameters->amp      ); 
    fprintf( info_file, "Sperm Number gamma= %e\n", parameters->Sperm_Number); 
    fprintf( info_file, "Dimensionless force B = %e\n", B); 

    tinercia = 1./parameters->gamma;
    Fmax = parameters->b_l*parameters->gamma/parameters->dt;

    fprintf( info_file, "tinercia = %e\n", tinercia); 
    if (tinercia < 10*parameters->dt){
      fprintf( info_file, "too big tinercia!\n");
    }
    fprintf( info_file, "Fmax = %e\n", Fmax); 
    if (parameters->Fx/parameters->mass > 0.1*parameters->b_l/parameters->dt){
      fprintf( info_file, "too big Fmax! \n"); 
    }

    if (parameters->oseen == 1){
      fprintf( info_file, "Oseen \n");
    }
    if (parameters->blake == 1){
      fprintf( info_file, "Blake \n");
    }
    fclose(info_file);
}

/* ########################################################################## */

void  simulate(     param_s  params,
                     vec_s  *r,
                     vec_s  *v)
/* ########################################################################## */
{
    vec_s *fc_1, *fc_2, *fc_3, *v_sv, *f_sv;
    vec_s *f_bend, *fc_1a, *fc_2a,*fc_3a;



    int    i, j, k, l,m, n_mon, step, n_ang, n_fil;

    float  *mags, *mags1;
    float  fx_tot, fy_tot, fz_tot, u, ke;
    float  rijx, rijy, rijz;
    float  work, fpar_ave, fper_ave;
    float  x_start, y_start, cycle_time, ang_tot;
    vec_s *rij, *rij1, *f,*fa, *tang,*v1, *r1, *rm, *rcm;

    FILE *conf = fopen("conf", "w");
   // FILE *flog_fa;
    /* Fab */ float freq, amplitude;

    /* Fab */ float fm_work, TOT_fm_work, AV_TOT_fm_work;
    /* Fab */ float ref_work, TOT_ref_work, AV_TOT_ref_work;
    /* Fab */ float velocity, TOT_velocity, AV_TOT_velocity, velocityx,velocityy, velocityz, velx, vely,velz;
    /* Fab */ int   count, TOT_count;
    /* Fab */ float Ext_force;
    /* Fab */ float Sperm_Number;
    /* Fab */ float x_old, y_old, z_old;
    /* Fab */ char filename[100];
    /* Fab */ FILE *file_ptr;

    /* MCL */ float  loc_dr, fatt, fatt2, b_l , y_amp,y_ampM,y_ampm;
    /* MCL */ float  dr0x, dr0y,drijx,drijy;
    /* MCL */ float dist;

    /* Fab */ vec_s constant_force;
    /* Fab */ constant_force.x = 0.0;
    /* Fab */ constant_force.y = 0.0;
    /* Fab */ constant_force.z = 0.0;

#ifdef LONG
    constant_force.x = params.Fx/params.mass;
#endif
#ifdef TRANS
    constant_force.z = params.Fz/params.mass;
    //constant_force.y = -params.gamma*params.b_l/params.dt/100.;
    printf(" constant_force.z = %f\n",constant_force.z);
#endif

    sprintf(filename,"%sGL_DR",params.dir_name);
   // flog_fa = fopen("fa.log","w");
       
    /* Fab */ freq = params.freq;
    /* Fab */ Sperm_Number = params.Sperm_Number;
    /* Fab */ amplitude = params.amplitude;
    /* Fab */ y_ampM = 0;
    /* Fab */ y_ampm = 0;
    /* Fab */ y_amp = 0;

    // printf("AMPLITUDE = %f\n",amplitude);
    // printf("cycle vx             vy              work       ref_work\n");    
    
    dist  = params.dist;
    n_fil = params.n_fil;
    n_mon = params.n_mon;   
    
    
/* Forze  a = di tutti  */    

    f_bend   = (vec_s *)malloc( n_mon*sizeof( vec_s ));

    fc_1a    = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));  
    fc_1     = (vec_s *)malloc( n_mon*sizeof( vec_s ));

    fc_2a    = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));
    fc_2     = (vec_s *)malloc( n_mon*sizeof( vec_s ));
 
    fc_3a    = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));
    fc_3     = (vec_s *)malloc( n_mon*sizeof( vec_s ));

    v_sv     = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));
    f_sv     = (vec_s *)malloc( n_mon*n_fil*sizeof( vec_s ));
    
    

    rij   = (vec_s *)malloc((n_mon + 1)*n_fil*sizeof( vec_s ));
    rij1  = (vec_s *)malloc((n_mon + 1)*sizeof( vec_s ));

    fa    = (vec_s *)malloc(n_mon*n_fil*sizeof( vec_s ));
    f     = (vec_s *)malloc(n_mon*sizeof( vec_s ));

    v1    = (vec_s *)malloc(n_mon*sizeof( vec_s ));
    r1    = (vec_s *)malloc(n_mon*sizeof( vec_s ));

    mags  = (float *)malloc(n_mon*n_fil*sizeof( float ));
    mags1 = (float *)malloc(n_mon*sizeof( float ));
    tang  = (vec_s *)malloc(n_mon*sizeof( vec_s ));
  
    rm    = (vec_s *)malloc(n_mon*sizeof( vec_s ));
    rcm = (vec_s *)malloc(n_fil*sizeof( vec_s ));


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


    //tang solo del primo
        tang[0].x       = -rij[1].x/mags[1];
	tang[0].y       = -rij[1].y/mags[1];
	tang[0].z       = -rij[1].z/mags[1];
	tang[n_mon-1].x = -rij[n_mon-1].x/mags[n_mon-1];
	tang[n_mon-1].y = -rij[n_mon-1].y/mags[n_mon-1];
	tang[n_mon-1].z = -rij[n_mon-1].z/mags[n_mon-1];
	for( i = 1; i < n_mon - 1; i++ )
          {
	    tang[i].x = -0.5*(rij[i].x/mags[i] + rij[i+1].x/mags[i+1]);
	    tang[i].y = -0.5*(rij[i].y/mags[i] + rij[i+1].y/mags[i+1]);
	    tang[i].z = -0.5*(rij[i].z/mags[i] + rij[i+1].z/mags[i+1]);
          }
    //endtang



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
      for( i = 0; i < n_mon; i++ ){
	rij1[i].x = rij[i+m*n_mon].x;
	rij1[i].y = rij[i+m*n_mon].y;
	rij1[i].z = rij[i+m*n_mon].z;


	mags1[i] = mags[i+m*n_mon];
	
	f[i].x = fa[i+m*n_mon].x;
	f[i].y = fa[i+m*n_mon].y;
	f[i].z = fa[i+m*n_mon].z;
      }
	
      u    = bending_forces( rij1, mags1, f, params );
     
      for( i = 0; i < n_mon; i++ ){
	fa[i+m*n_mon].x = f[i].x;
	fa[i+m*n_mon].y = f[i].y;
	fa[i+m*n_mon].z = f[i].z;
	
      }
    }

#ifdef CONSTANT_FORCE
    /* Fab */
    /* Add a constant Force all over the flagella */
    /* */for( i = 0; i < n_mon; i++ )
    /* */{
    /* */   fa[i].x += constant_force.x;
    /* */   fa[i].y += constant_force.y;
    /* */   fa[i].z += constant_force.z;
    /* */ }  
    /* end Fab */
#endif
   
#ifdef HEAD
    /* Fab */
    /* Add an oscillating force on the head */
    /* There is a "k" which here is = 0 (is the time....)*/
    /* */   fa[0].y += amplitude*cos(freq*0.0*params.dt);
    /* */   fa[n_mon].y += amplitude*cos(freq*0.0*params.dt);
    /* end Fab */
#endif


#ifdef TORQUE_M
	     fatt       =  amplitude*sin(freq*0.0*params.dt);     
	  
	     
	     dr0x       =   (r[1].x + r[0].x)/2.0;
	     dr0y       =   (r[1].y + r[0].y)/2.0;
	     drijx      =  -(r[0].x - dr0x);
	     drijy      =  +(r[0].y - dr0y);
	     
	     fa[0].x += drijy*fatt;
	     fa[0].y += drijx*fatt;
	     drijx      =  -(r[1].x - dr0x);
	     drijy      =  +(r[1].y - dr0y);
	     fa[1].x   += drijy*fatt;
	     fa[1].y   += drijx*fatt;
#endif

	     /* Fab */
	     TOT_ref_work = 0.0;
	     TOT_fm_work  = 0.0;
	     TOT_velocity = 0.0;
	     TOT_count    = 0;
	     
 

	     /*  inizio ciclo temporale */

////////////////////////////////////////////////////////////////////////////////////////////
//COMENÇA EL CICLE TEMPORAL!!!!!!!!!!!!!!!
////////////////////////////////////////////////////////////////////////////////////////////

   for( j = 0; j < params.n_steps; j++ )
    {
       x_start  = r[n_mon-1].x;
       y_start  = r[n_mon-1].y;
       work     = 0.0;
       n_ang    = 0;
       ang_tot  = 0.0;       
     
       /* Fab */
       fm_work  = 0.0;
       velocity = 0.0;
       velocityx = 0.0;
       velocityy = 0.0;
       velocityz = 0.0;
       count    = 0;

       for( i = 0; i < n_mon; i++ ){
	  rm[i].x=0.0;
	  rm[i].y=0.0;
	  rm[i].z=0.0;}
      
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
	    for( i = 0; i < n_mon; i++ ){
	      rij1[i].x = rij[i+m*n_mon].x;
	      rij1[i].y = rij[i+m*n_mon].y;
	      rij1[i].z = rij[i+m*n_mon].z;
	 
	   
	      mags1[i] = mags[i+m*n_mon];
	    
	      v1[i].x = v[i+m*n_mon].x;
	      v1[i].y = v[i+m*n_mon].y;
	      v1[i].z = v[i+m*n_mon].z;
	
	    }
	     

	   
	     constrain_velocities( rij1, mags1, v1,fc_1, params );
	
	    for( i = 0; i < n_mon; i++ ){
	      v[i+m*n_mon].x = v1[i].x;
	      v[i+m*n_mon].y = v1[i].y;
	      v[i+m*n_mon].z = v1[i].z;
	     
	

	      fc_1a[i+m*n_mon].x   = fc_1[i].x;
	      fc_1a[i+m*n_mon].y   = fc_1[i].y;
	      fc_1a[i+m*n_mon].z   = fc_1[i].z;

	    }
	  }
	   
	
          for( i = 0; i < n_mon*n_fil; i++ )
          {
              v_sv[i].x = v[i].x;
              v_sv[i].y = v[i].y;
              v_sv[i].z = v[i].z;
          }       
	  //fprintf(flog_fa,"  b: %le",fa[15].y);
          ke = verlet_pt1( r, v, fa, tang, params );


	
	  for(m=0; m<n_fil; m++){
	    for( i = 0; i < n_mon; i++ ){
	      rij1[i].x = rij[i+m*n_mon].x;
	      rij1[i].y = rij[i+m*n_mon].y;
	      rij1[i].z = rij[i+m*n_mon].z;
	      
	      mags1[i] = mags[i+m*n_mon];

	      v1[i].x = v[i+m*n_mon].x;
	      v1[i].y = v[i+m*n_mon].y;
	      v1[i].z = v[i+m*n_mon].z;

	
	      r1[i].x = r[i+m*n_mon].x;
	      r1[i].y = r[i+m*n_mon].y;
	      r1[i].z = r[i+m*n_mon].z;
	
	    }    
	    constrain_positions( r1, rij1, v1, fc_2, params );
	    
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


	

	  //tamg boh
          tang[0].x       = -rij[1].x/mags[1];
          tang[0].y       = -rij[1].y/mags[1];
          tang[0].z       = -rij[1].z/mags[1];
          tang[n_mon-1].x = -rij[n_mon-1].x/mags[n_mon-1];
          tang[n_mon-1].y = -rij[n_mon-1].y/mags[n_mon-1];
          tang[n_mon-1].z = -rij[n_mon-1].z/mags[n_mon-1];
          for( i = 1; i < n_mon - 1; i++ )
          {
             tang[i].x = -0.5*(rij[i].x/mags[i] + rij[i+1].x/mags[i+1]);
             tang[i].y = -0.5*(rij[i].y/mags[i] + rij[i+1].y/mags[i+1]);
             tang[i].z = -0.5*(rij[i].z/mags[i] + rij[i+1].z/mags[i+1]);
          }
	  //tang boh

          for( i = 0; i < n_mon*n_fil; i++ )
          {
	
	    fa[i].x = 0.0;
	    fa[i].y = 0.0;
	    fa[i].z = 0.0;
	  }

	 
	  
#ifdef HEAD
	  /* Fab */
	  /* Add an oscillating force on the head */
          if(y_ampm > r[0].y) y_ampm = r[0].y;
	  if(y_ampM < r[0].y) y_ampM = r[0].y;
//	  printf("%f %f \n", y_ampm,y_ampM); 

	  Ext_force = amplitude*cos(freq*k*params.dt);
	  fa[0].y += Ext_force;
	  fa[n_mon].y += Ext_force;
	  //   f[0].y += (Ext_force=amplitude*cos(freq*k*params.dt)); 
	  /* end Fab */
#endif  

#ifdef CONSTANT_FORCE
          /* Fab */
          /* Add a constant Force all over the flagella */
	  /*  for( i = 0; i < n_mon*n_fil; i++ ) */
          for( i = 0; i < n_mon; i++ )
          {
             fa[i].x += constant_force.x;
             fa[i].y += constant_force.y;
             fa[i].z += constant_force.z;}
#endif


	  //questi sono da trasf in pot armonici (trap)
	  //fprintf(flog_fa,"  c: %le",fa[15].y);
#ifdef TORQUE_M
	     fatt       =  amplitude*sin(freq*k*params.dt);     
	  
	     
	     dr0x       =   (r[1].x + r[0].x)/2.0;
	     dr0y       =   (r[1].y + r[0].y)/2.0;
	     drijx      =  -(r[0].x - dr0x);
	     drijy      =  +(r[0].y - dr0y);
	     
	     fa[0].x += drijy*fatt;
	     fa[0].y += drijx*fatt;
	     drijx      =  -(r[1].x - dr0x);
	     drijy      =  +(r[1].y - dr0y);
	     fa[1].x   += drijy*fatt;
	     fa[1].y   += drijx*fatt;
#endif




#ifdef TORQUE_PINNED
	     b_l = params.b_l;
	     fatt       =  amplitude*sin(freq*k*params.dt);     
	     fatt2      =  amplitude*freq*cos(freq*k*params.dt);
	
	     
	     r[0].y   = 0;
	     v[0].y   = 0;
	     
       	     
	     r[1].y   = b_l*fatt;
	     v[1].y   = b_l*fatt2;	     
#endif

#ifdef TORQUE_PINNED_FIX
	     b_l = params.b_l;
	     fatt       =  amplitude*sin(freq*k*params.dt);     
	     fatt2      =  amplitude*freq*cos(freq*k*params.dt);
	
	     
	     r[0].x   = 0;
	     r[0].y   = 0;
	     v[0].x   = 0;
	     v[0].y   = 0;
	     
       	     r[1].x   = b_l - b_l*fatt;
	     r[1].y   = b_l*fatt;
	     v[1].x   = b_l*fatt2;
	     v[1].y   = - b_l*fatt2;
	     

#endif

	   
	    // fprintf(flog_fa,"  d1: %le",fa[15].y);
	     //work += boundary_conditions( r, v, fa, tang, params, step, j);
	     //fprintf(flog_fa,"  d2: %le",fa[15].y);
#ifdef WAVE          
          if(y_ampM < r[0].y) y_ampM = fabs(r[0].y);

	  y_ampm = -y_ampM;
	 
#endif 
	  
	  for(m=0; m<n_fil; m++){
	    for( i = 0; i < n_mon; i++ ){

	      
	      f_bend[i].x = 0.0;
	      f_bend[i].y = 0.0;
	      f_bend[i].z = 0.0;
	

	      rij1[i].x = rij[i+m*n_mon].x;
	      rij1[i].y = rij[i+m*n_mon].y;
	      rij1[i].z = rij[i+m*n_mon].z;


	      mags1[i] = mags[i+m*n_mon];

	
	    }
     
	    u  = bending_forces( rij1, mags1, f_bend, params );
     

	    for( i = 0; i < n_mon; i++ ){
	      fa[i+m*n_mon].x += f_bend[i].x;
	      fa[i+m*n_mon].y += f_bend[i].y;
	      fa[i+m*n_mon].z += f_bend[i].z;	
      
	    }
	  }
	 // fprintf(flog_fa,"  e: %le",fa[15].y);
	 
/*    for( i = 0; i < n_mon*n_fil; i++ ) */
/*            { */
/*                f_sv[i].x = fa[i].x; */
/*                f_sv[i].y = fa[i].y; */
/*                f_sv[i].z = fa[i].z; */
/*     	      } */

	  for(m=0; m<n_fil; m++){
	    for( i = 0; i < n_mon; i++ ){
	      rij1[i].x = rij[i+m*n_mon].x;
	      rij1[i].y = rij[i+m*n_mon].y;
	      rij1[i].z = rij[i+m*n_mon].z;
	 
	   
	      mags1[i] = mags[i+m*n_mon];
	    
	      f[i].x = fa[i+m*n_mon].x;
	      f[i].y = fa[i+m*n_mon].y;
	      f[i].z = fa[i+m*n_mon].z;
	
	    }
	     

	   
	      constrain_forces(rij1, mags1, f, fc_3, params );
	
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
	    
	    f_sv[i].x = fa[i].x;
	    f_sv[i].y = fa[i].y;
	    f_sv[i].z = fa[i].z;
	      
          }       
	 // fprintf(flog_fa,"  f: %le\n",fa[15].y);
	  verlet_pt2( r, v, fa ,tang, params );

      /* Fab */
	  if((k%10==0)&&(j>2))
	    {
	      fm_work += fabs(params.mass*(f[0].y+(fc_1[0].y+fc_2[0].y)/2)
			      *(r[0].y - y_old));   
	        
	      count++;
	    }
	  /* End Fab */

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
	

	  for( i = 0; i < n_mon; i++ ){
	  rm[i].x += r[i].x;
	  rm[i].y += r[i].y;
	  rm[i].z += r[i].z;}


	  

    }
    
/*  fine ciclo temporale */
/*	  for( m = 0; m < n_fil; m++ ){
	 for(i=0; i<n_mon; i++)
	   rcm[m].z +=  r[i+m*n_mon].z;
	 rcm[m].z /= n_mon;}
*/
	  /*
       printf("old\trz0=%e\tzN=%e\trcmz0=%e\trcmz1=%e\t", r[0].z, r[n_mon].z, rcm[0].z, rcm[1].z);
	  
	       for( m = 0; m < n_fil; m++ ){
	 for(i=0; i<n_mon; i++)
	   r[i+m*n_mon].z -= rcm[m].z - m*dist;}
	       
       printf("new\trz0=%e\trzN=%e\n", r[0].z, r[n_mon].z);
       */

       /* Fab */
       y_amp = fabs(y_ampM - y_ampm)/2;
       
      
       
     /* Fab */ /* Accumulates average over cycles */
       velocityx /= params.n_cycle;
       velocityy /= params.n_cycle;
       velocityz /= params.n_cycle;
       velocity = sqrt(velocityx*velocityx+velocityy*velocityy+velocityz*velocityz);
	  
       	 //if(j==params.n_steps-1)   
	 printf("%i\t%le\t%le\t%le\t%le\t%le\n", j, velocity, velocityx,velocityy, velocityz, y_amp);
       
       velocityx = 0.0;
       velocityy = 0.0;
       velocityz = 0.0;
	

//       velocity = (r[n_mon-1].z - z_start)/cycle_time;
       
       for( i = 0; i < n_mon; i++ ){
	 rm[i].x /= params.n_cycle;
	 rm[i].y /= params.n_cycle;
	 rm[i].z /= params.n_cycle;}
       
       
       output_data( rm, v, fa, params, step, &fper_ave, &fpar_ave );
       fflush( stdout );
	  


    
       ref_work = params.gamma * params.mass * n_mon;
       ref_work *= velocity*velocity;

       if(j>2){ 
	 TOT_fm_work  += fm_work/count;
	 TOT_velocity += velocity;
	 TOT_ref_work += ref_work;
	 
	 TOT_count++;
       }
       // printf(" %12.5g %12.5g\n", fm_work/count, ref_work);
       /* End Fab */

   }
 //  fclose(flog_fa);
/*  fine tutti cicli  */

     /* FM: Average over all cycles */
    AV_TOT_fm_work  = TOT_fm_work/TOT_count;
    AV_TOT_ref_work = TOT_ref_work/TOT_count;
    AV_TOT_velocity = TOT_velocity/TOT_count;

  
    sprintf(filename,"%sAV_VELOCITY",params.dir_name);
    file_ptr = fopen(filename,"w");
    fprintf(file_ptr,"%.16g %.16g %.16g %.16g\n"
	    ,Sperm_Number, (AV_TOT_velocity/(y_amp*y_amp)),
	    AV_TOT_velocity, y_amp); 
    fclose(file_ptr);

    sprintf(filename,"%sAV_WORK",params.dir_name);
    file_ptr = fopen(filename,"w");
    fprintf(file_ptr,"%.16g %.16g %.16g %.16g\n",Sperm_Number,
	    AV_TOT_fm_work, AV_TOT_ref_work,
	    AV_TOT_ref_work/AV_TOT_fm_work);   
    fclose(file_ptr);

    sprintf(filename,"%sAV_EFFICIENCY",params.dir_name);
    file_ptr = fopen(filename,"w");
    fprintf(file_ptr,"%.16g %.16g\n",Sperm_Number,
	    AV_TOT_ref_work/AV_TOT_fm_work);  
    fclose(file_ptr);

#ifdef TRANS
    sprintf(filename,"%sGL_DR",params.dir_name);
    file_ptr = fopen(filename,"w");
    for( l = 0; l < n_mon; l++ )
      {	 
	loc_dr  = (float)((constant_force.z * params.mass)/ (v[l].z));
	fprintf(file_ptr,"%.16g %.16g\n", 
		(float)l/(float)(n_mon-1),(n_mon)* loc_dr);  
      } 
    fclose(file_ptr);
#endif
    
#ifdef LONG
    sprintf(filename,"%sGL_DR",params.dir_name);
    file_ptr = fopen(filename,"w");
    for( l = 0; l < n_mon; l++ )
      {	 
	loc_dr  = (float)((constant_force.x * params.mass)/ (v[l].x));
	fprintf(file_ptr,"%.16g %.16g\n",  
		(float)l/(float)(n_mon-1), n_mon* loc_dr);  
      } 
    fclose(file_ptr);
#endif
   

}



