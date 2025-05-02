//Function for remeshing:
#include "M4.h"
#include <math.h>
//Function code:
double M4(double x_tilda, double x, double h )
{
double u= fabs(x_tilda - x)/(h);
double value = 0.0;
if( u>= 0.0 && u <=1.0)
{
value = 1- 2.5*u*u + 1.5*u*u*u;
return value;
}
else if( u>=1.0 && u<=2.0)
{
value = 1.0/2.0*(2.0-u)*(2.0-u)*(1.0-u);
return value;
}
else
{
return 0.0;
}
}