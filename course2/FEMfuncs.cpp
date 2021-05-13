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

real FEM::theta(knot point, real t)
{
	return -point.x * t;
}

real FEM::f(knot point, real t)
{
	switch (un)
	{
	case 1: return 0;
	case 2: return -6. * point.y;
	case 3: return point.y;
	case 4: return 3 * t * t + 6 * t;
	case 5: return 2 * t + 2;
	case 6: return 2 * t - 2;
	case 7: return point.x * point.y;
	case 8: return 0;
	case 9: return 1;
	case 10: return 2 * t + 2;
	case 11: return 6 * t * t + 3 * t;

	default:
		return 0;
	}
}

real FEM::ubeta(knot point, real t)
{
	return ug(point, t) - point.x * t;
}

real FEM::ug(knot point, real t)
{
	switch (un)
	{
	case 1: return point.x * point.y + 1;
	case 2: return point.y * point.y * point.y; 
	case 3: return point.y * t;
	case 4: return point.x + point.y + t * t * t;
	case 5: return point.x + point.y + t * t;
	case 6: return point.x * point.x + point.y * point.y + t * t;
	case 7: return point.x * point.y * t;
	case 8: return 1;
	case 9: return t;
	case 10: return t * t;
	case 11: return t * t * t;

	default:
		return 0;
	}
}

