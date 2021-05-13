#include "FEM.h"
#include <math.h>

real FEM::lambda(knot point, real t)
{
	switch (un)
	{
	

	default:
		return 1;
	}
}

real FEM::f(knot point, real t)
{
	switch (un)
	{
	case 1: return -2. * t + point.y * point.y;
	case 2: return -6. * point.y;
	case 3: return point.y;
	case 4: return 3 * t * t + 6 * t;
	case 5: return 2 * t + 2;
	case 6: return 2 * t - 2;
	case 7: return point.x * point.y;

	default:
		return 0;
	}
}

real FEM::ug(knot point, real t)
{
	switch (un)
	{
	case 1: return point.y * point.y * t + 1;
	case 2: return point.y * point.y * point.y; 
	case 3: return point.y * t;
	case 4: return point.x + point.y + t * t * t;
	case 5: return point.x + point.y + t * t;
	case 6: return point.x * point.x + point.y * point.y + t * t;
	case 7: return point.x * point.y * t;

	default:
		return 0;
	}
}

