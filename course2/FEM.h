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
      bool bound_num;
      int bound_param;
   };
   list<bound*> bounds;

   struct Matrix {
      real *l, *u, *d;
      int *ig, *jg;
   };

   knot* knots;

   //real lambda(real knot[2], real t);
   real f(knot knot_, real t);
   real gamma = 1, theta = 1, khi = 1, lambda = 1;
   int num_of_knots, num_of_FE, un, nt;

   real localM[4][4]; // 4*4
   real localG[4][4];
   real localA[4][4];

   void MakeSparseFormat();
   void AddElement(Matrix *A, int knot_num[4], int i, int j, real elem);
   void AddLocal(Matrix* A, int knot_num[4], real localA[4][4], real coeff);
   void MakeBounds(int s);
   void CreateA(int s);
   void CreateM(rectangle *rect); // можно посылать параметр, на кот. умножаетс€ матрица
   void CreateG(rectangle *rect); // л€мбду надо при формировании
   void Createb(int s);
   void Output();
   real diffNorm(real *q);
   real norm(real* b);
   real ug(knot knot_, real t);
   //real **A, **M, **G;
   Matrix *A, *M, *G;
   real *b, *qk, *temp, *t;
   real **q;  // q[0] & q[1] - начальные
   void equalize(real *x, real *y);
   real scalar(real* v, real* u);
   void MulAb(real* v, Matrix* A, real* b);
   real* z, *r, *p, *ff, *x;
   void SolveSLAE(int s);

   public:
   FEM();
   void Solve();
   void ShowError(int s, ofstream &out);    
};