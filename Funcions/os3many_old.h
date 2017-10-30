typedef struct {
int   n_steps;
int   n_mon;
int   n_cycle;
int   n_fil;/*   Mar */
int   blake;
int   oseen;
float mass;
float gamma;
float hydr_radius;
float gamma_head;
float dt;
float k;
float phi_0;
float b_l;
float tolerance;
float amp;
float freq;
float w_length;
int   n_left;
int   n_right;
int   out_freq;
char  infile_name[MAX_BUFF];
char  outfile_name[MAX_BUFF];
char  dir_name[MAX_BUFF];
float amplitude; /* Fab */
float Sperm_Number; /* Fab */
float epsD0;
float epsBferro0; 
float dist; /*  Mar */
float Fx;
float Fz;
float Bx;
float Bz;
float chi;
float m2;
} param_s;


void   configure_sys( int,
                      char    *[],
                      param_s *,
                      vec_s   **,
                      vec_s   **);

void   simulate( param_s,
                 vec_s *,
                 vec_s *  );

/* Dipolar force routine
Compute dipolar force between
paramagnetic and ferromagnetic 
beads at the edges of the filament

Return: f_dip forces
*/

void  dipolar_forces( vec_s,
		       vec_s,
		       vec_s *,
		       param_s,
		       vec_s);

/*  Torque Bfield routine
Compute couple of forces caused
by the external field on the ferromagnet

Return: f_Btorque forces
*/

void  Torque_Bfield( vec_s ,
           vec_s  ,
           vec_s  ,
           vec_s  *,
           param_s,
           vec_s);

/* Bending force routine
Compute bending forces. 
It calls angle_force
Return:   u (???),
          f (forces)
          others?
*/
float  bending_forces( vec_s *,
                       float *,
                       vec_s *,
                       param_s  );
/**************************
/* Bending force routine
Compute bending forces. 
It calls angle_force
Return:   u (???),
          f (forces)
          others?
*/
float  bending_forces_circ( vec_s *,
                       float *,
                       vec_s *,
                       param_s  );
/**************************
Computes angular forces,
diedral????
The forces are combinations of distances......????
*/
float  angle_force( vec_s,
                    vec_s,
                    float,
                    float,
                    float,
                    float,
                    vec_s * );

void    constrain_positions(  vec_s *,
                              vec_s *,
                              vec_s *,
                              vec_s *,
                              param_s );

void   constrain_velocities(  vec_s *, 
                              float *,
                              vec_s *,
                              vec_s *,			      
                              param_s );

void   constrain_forces(      vec_s *, 
                              float *,
                              vec_s *,
                              vec_s *,
                              param_s );

void    constrain_positions_circ(  vec_s *,
                              vec_s *,
                              vec_s *,
                              vec_s *,
                              param_s );

void   constrain_velocities_circ(  vec_s *, 
                              float *,
                              vec_s *,
                              vec_s *,			      
                              param_s );

void   constrain_forces_circ(      vec_s *, 
                              float *,
                              vec_s *,
                              vec_s *,
                              param_s );

float   verlet_pt1( vec_s *,
                    vec_s *,
                    vec_s *,
                    vec_s *,
                    param_s );



void   verlet_pt2( vec_s *,
                   vec_s *,
                   vec_s *,
                   vec_s *,
                   param_s );

float   verlet_pt1_gamma_var( vec_s *,
                    vec_s *,
                    vec_s *,
                    vec_s *,
                    param_s );



void   verlet_pt2_gamma_var( vec_s *,
                   vec_s *,
                   vec_s *,
                   vec_s *,
                   param_s );
void oseen( vec_s *,
	    vec_s *,
	    vec_s *,
	    param_s);

void blake( vec_s *,
	    vec_s *,
	    vec_s *,
	    param_s);

void  output_data( vec_s *,
                   vec_s *,
                   vec_s *,
                   param_s,
                   int,
                   float *,
                   float *);


float  boundary_conditions( vec_s *,
                            vec_s *,
                            vec_s *,
                            vec_s *,
                            param_s,
                            int,
                            int        );

/* Fab */
void print_force(vec_s *, param_s, int, char *);
