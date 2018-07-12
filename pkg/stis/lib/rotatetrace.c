#include <stdio.h>
#include <stdlib.h>
#include<math.h>

# define PI   3.14159265358979323846

double rotatetrace(double expstart, double mjd, double degperyr, double *a2displ, int nelem) {
  
  int i;
  double temp;
  double angle = degperyr * (expstart - mjd) / 365.25;
  double radians_coef = PI / 180;
 
  temp = tan(angle * radians_coef);
  for (i = 0; i < nelem; i++){
   a2displ[i] -= (i - nelem/2) * temp;
  }    


  return angle;
}
