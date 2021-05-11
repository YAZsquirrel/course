#include "FEM.h"
#include <math.h>

//real FEM::lambda(real knot[2], real t)
//{
//	return real();
//}

real FEM::f(knot point, real t)
{
	switch (un)
	{
	case 1: return 2. + point.y * point.y;
	case 2: return -6. * point.y + point.y * point.y * point.y;
	case 3: return point.y;
	case 4: return 0;
	case 5: return point.x * point.x + point.y * point.y + point.x * point.y + 1 - 4;

	default:
		return 0;
	}
}

real FEM::ug(knot point, real t)
{
	switch (un)
	{
	case 1: return point.y * point.y;
	case 2: return point.y * point.y * point.y; 
	case 3: return point.y;
	case 4: return 0;
	case 5: return point.x * point.x  + point.y * point.y  + point.x * point.y + 1;

	default:
		return 0;
	}
}
