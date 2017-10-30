#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//#include "nr.h"
#include "nrutil.h"
#include "nrutil.c"
#include "tridag.c"
//#define float double

void main()
{
  static float *a, *b, *c, *rhs, *ans;
  a      = (float *)malloc( (6)*sizeof( float )); 
  b      = (float *)malloc( (6)*sizeof( float )); 
  c      = (float *)malloc( (6)*sizeof( float )); 
  ans    = (float *)malloc( (6)*sizeof( float )); 
  rhs    = (float *)malloc( (6)*sizeof( float ));
/*  a[0] = 0.0;
  b[0] = 1.0;
  c[0] = 2.0;
  a[1] = 3.0;
  b[1] = 4.0;
  c[1] = 5.0;
  a[2] = 6.0;
  b[2] = 7.0;
  c[2] = 0.0;
  
  rhs[0] = 5.0;
  rhs[1] = 26.0;
  rhs[2] = 33.0;
*/

  a[1] = 0.0;
  b[1] = 1.0;
  c[1] = 2.0;
  
  a[2] = 3.0;
  b[2] = 4.0;
  c[2] = 5.0;
  
  a[3] = 6.0;
  b[3] = 7.0;
  c[3] = 0.0;
 
  
//  rhs[0] = 0.0;
  rhs[1] = 5.0;
  rhs[2] = 26.0;
  rhs[3] = 33.0;


  
  tridag(a,b,c,rhs, ans,3);
  printf("ans0 = %f, ans1 = %f, ans2 = %f, ans3 = %f  \n",ans[0], ans[1], ans[2], ans[3]);
 
}