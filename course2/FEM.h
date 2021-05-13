#pragma once
#include <iostream>
#include <list>
typedef double real;
using namespace std;

class FEM
{
   private:
   struct knot
   {
      real x, y;
   };

   struct rectangle {
      real hy = 0, hx = 0; 
      int knots_num[4]{};
      knot knots[4]{};
   };
   //list<rectangle*> rects;
   rectangle** rects;

   struct bound {
      int knot_nums[2];
      knot knots[2];
      int bound_num;
      int bound_param;
   };
   list<bound*> bounds;

   struct Matrix {
      real *l, *u, *d;
      int *ig, *jg;
   };

   knot* knots;

   real lambda(knot knot_, real t);
   real theta(knot knot_, real t);
   real f(knot knot_, real t);
   real ubeta(knot knot_, real t);
   real sigma = 1, khi = 1;
   int num_of_knots, num_of_FE, un, nt, num_of_not_zero;

   real localM[4][4]; // 4*4
   real localG[4][4];
   real localA[4][4];

   void MakeSparseFormat();
   void AddElement(Matrix *A, int i, int j, real elem);
   void AddLocal(Matrix* A, int knot_num[4], real localA[4][4], real coeff);
   void MakeBounds(int s);
   void CreateA(int s);
   void CreateM(rectangle *rect); // можно посылать параметр, на кот. умножается матрица
   void CreateG(rectangle *rect); // лямбду надо при формировании
   void Createb(int s);
   void equalize(real *x, real *y);
   void SolveSLAE(int s);
   void MulAb(real* v, Matrix* A, real* b);
   void nullify(real *x);
   real scalar(real* v, real* u);
   real ug(knot knot_, real t);

   real* z, *r, *p, *ff;
   Matrix *A, *M, *G;
   real *b, *b0, *temp, *t;
   real *q0, *q1, *q;  // q0 & q1 - начальные, а также q0 = q^k-2, q1 = q^k-1

   public:
   FEM();
   void Solve();
   void Output(int s, ofstream &out);    
};