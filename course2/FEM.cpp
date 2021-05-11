#include "FEM.h"
#include <fstream>
#include <iostream>

FEM::FEM()
{
#pragma region input

   ifstream fknots("Knots.txt");
   ifstream frects("Rectangles.txt");
   ifstream fbounds("Bounds.txt");
   ifstream fparams("Params.txt");

   fknots >> num_of_knots;
   knots = new knot[num_of_knots];
   for (int i = 0; i < num_of_knots; i++)
      fknots >> knots[i].x >> knots[i].y;
   fknots.close();

   rectangle *rect;
   frects >> num_of_FE;
   rects = new rectangle*[num_of_FE];
   int num;
   for (int i = 0; i < num_of_FE; i++)
   {
      rect = new rectangle();
      for (int k = 0; k < 4; k++)
      {
         frects >> num;
         rect->knots_num[k] = num - 1;
         rect->knots[k] = knots[num - 1];
      }
      rect->hy = rect->knots[3].y - rect->knots[0].y;
      rect->hx = rect->knots[3].x - rect->knots[0].x;
      rects[i] = rect; 
      //rects.push_back(rect);
   }
   frects.close();
   
   int cond_num;
   fbounds >> cond_num;
   for (int i = 0; i < cond_num; i++)
   {
      bound* cond = new bound;
      fbounds >> num;	//номер условия
      cond->bound_num = num == 1;
      if (!cond->bound_num)
        fbounds >> cond->bound_param;
      int number;
      for (int i = 0; i < 2; i++)
      {
         fbounds >> number;									// 
         cond->knot_nums[i] = number - 1;				// номера узлов ребра, на котором краевое...
         cond->knots[i] = knots[number - 1];	// 
      }
      //fcond >> number;										//
      bounds.push_back(cond);
   }
   fbounds.close();

   fparams >> lambda >> gamma >> khi >> un >> nt;
   t = new real[nt];
   for (int i = 0; i < nt; i++)
      fparams >> t[i];
   fparams.close();
#pragma endregion
   
   MakeSparseFormat();
   q = new real*[nt];
   for (int i = 0; i < nt; i++)
      q[i] = new real[num_of_knots]{};
   b = new real[num_of_knots]{};
   temp = new real[num_of_knots]{};
   ff = new real[num_of_knots] {};
   z = new real[num_of_knots] {};
   r = new real[num_of_knots] {};
   p = new real[num_of_knots] {};
   //qk = new real[num_of_knots]{};
}

void FEM::MakeSparseFormat() 
{
   const int N = 4; 
   int *list1, *list2;
   int *listbeg = new int[num_of_knots];

   for (int i = 0; i < num_of_knots; i++)
      listbeg[i] = -1;

   list1 = new int[num_of_knots * num_of_knots]{};
   list2 = new int[num_of_knots * num_of_knots]{};
   int listsize = -1, iaddr, ind1, ind2, k;

   for (int iel = 0; iel < num_of_FE; iel++) // проход по всем кэ
   {
      for (int i = 0; i < 4; i++) // по узлам кэ
      {
         k = rects[iel]->knots_num[i]; //
         for (int j = i + 1; j < 4; j++) // 
         {
            ind1 = k;
            ind2 = rects[iel]->knots_num[j];  //
            if (ind2 < ind1) //занесение связи большого номера с меньшим, те ind2 с ind1
            {
               ind1 = ind2;
               ind2 = k;
            }
            iaddr = listbeg[ind2];
            if (iaddr == -1) // если список пуст
            {
               listsize++;
               listbeg[ind2] = listsize;
               list1[listsize] = ind1;
               list2[listsize] = -1;
            }
            else // ищем в списке ind1
            {
               while (list1[iaddr] < ind1 && list2[iaddr] >= 0)
                  iaddr = list2[iaddr];
               if (list1[iaddr] > ind1)  // если не нашли и встретили эл с большим номером 
               {                         // добавляем перед ним, чтоб список был упорядоченным
                  listsize++;
                  list1[listsize] = list1[iaddr];
                  list2[listsize] = list2[iaddr];
                  list1[iaddr] = ind1;
                  list2[iaddr] = listsize;
               }
               else if (list1[iaddr] < ind1) // если не нашли, то пихаем в конец
               {
                  listsize++;
                  list2[iaddr] = listsize;
                  list1[listsize] = ind1;
                  list2[listsize] = -1;
               }
            }
         }
      }
   }

   A = new Matrix;
   G = new Matrix;
   M = new Matrix;
   G->ig = M->ig = A->ig = new int[num_of_knots + 1]{};
   G->jg = M->jg = A->jg = new int[listsize + 1]{};  // +1???

   for (int i = 0; i < num_of_knots; i++)
   {
      A->ig[i + 1] = A->ig[i];
      
      for (iaddr = listbeg[i]; iaddr != -1; )
      {
         A->jg[A->ig[i + 1]] = list1[iaddr];
         A->ig[i + 1]++;
         iaddr = list2[iaddr];
      }
   }

   delete[] listbeg;
   delete[] list1;
   delete[] list2;

   cout << "jg: ";
   for (int i = 0; i < listsize + 1; i++)
      cout << A->jg[i] << " ";
   cout << "\nig: ";
   for (int i = 0; i < num_of_knots + 1; i++)
      cout << A->ig[i] << " ";
   cout << "\n";

   A->l = new real[listsize + 1]{};
   M->l = new real[listsize + 1]{};
   G->l = new real[listsize + 1]{};

   A->u = new real[listsize + 1]{};
   M->u = new real[listsize + 1]{};
   G->u = new real[listsize + 1]{};

   A->d = new real[num_of_knots]{};
   M->d = new real[num_of_knots]{};
   G->d = new real[num_of_knots]{};

}

void FEM::AddElement(Matrix* A, int knot_num[4], int i, int j, real elem)
{
   if(i == j)
      A->d[knot_num[i]] += elem;
   else if (i < j)
   {
      for (int i = A->ig[j]; i < A->ig[j + 1]; i++)
         if (A->jg[i] == i) break;
      
      A->u[i] += elem; // i-1?
   }
   else
   {
      for (int i = A->ig[j]; i < A->ig[j + 1]; i++)
         if (A->jg[i] == j) break;

      A->l[i] += elem; // i-1??
   }

}

void FEM::AddLocal(int knot_num[4], real localA[4][4])
{
   int ibeg, iend, ind;
   for (int i = 0; i < 4; i++)
      A->d[knot_num[i]] += localA[i][i];
   for (int i = 0; i < 4; i++)
   {
      ibeg = A->ig[knot_num[i]];
      for (int j = 0; j < i; j++) // i - 1?
      {
         iend = A->ig[knot_num[i] + 1];  // -1 ?
         while (A->jg[ibeg] != knot_num[j])
         {
            ind = (ibeg + iend)/2;
            if (A->jg[ind] < knot_num[j])
               ibeg = ind + 1;
            else
               iend = ind;
         }
         A->l[ibeg] += localA[i][j];
         A->u[ibeg] += localA[j][i];
         ibeg++;
      }
   }
   
}

void FEM::equalize(real* x, real* y)
{
   for (int i = 0; i < num_of_knots; i++)
      x[i] = y[i];
}

void FEM::Solve()
{
   CreateA(0);
   Createb(0);
   MakeBounds(0);
   SolveSLAE(0);
   ShowError(0);
}

void FEM::ShowError(int s)
{
   for (int i = 0; i < num_of_knots; i++)
      cout << q[s][i] - ug(knots[i], s) << '\n';
}

void FEM::MakeBounds(int s)
{
   for (list<bound*>::iterator iter = bounds.begin(); iter != bounds.end(); iter++)
   {
      bound* cond = *iter;

      if (cond->bound_num)
      {
         for (int i = 0; i < 2; i++)
         {
            A->d[cond->knot_nums[i]] = 1;
            for (int j = A->ig[cond->knot_nums[i]] ; j < A->ig[cond->knot_nums[i] + 1]; j++)
               A->l[j] = 0.;
            for (int j = 0; j < A->ig[num_of_knots]; j++)
               if (A->jg[j] == cond->knot_nums[i])
                  A->u[j] = 0.;
           
            b[cond->knot_nums[i]] = ug(cond->knots[i], s);
         }
      }
      else
      {

      }
   }

}

void FEM::CreateA(int s) // 1/dt^2 Mx + 1/2dt Mo + 1/2 G
{
   
   rectangle *rect;
   for (int i = 0; i < num_of_FE; i++)
   {
      rect = rects[i];
      CreateG(rect);
      CreateM(rect);
      AddLocal(rect->knots_num, localM);
      AddLocal(rect->knots_num, localG);
   }

}

void FEM::CreateM(rectangle *rect)
{
   localM[0][0] = localM[1][1] = localM[2][2] = localM[3][3] = gamma / 9. * (rect->hx) * (rect->hy); //диагональ
   localM[0][1] = localM[1][0] = localM[0][2] = localM[2][0] =
   localM[2][3] = localM[3][2] = localM[1][3] = localM[3][1] = gamma / 18. * (rect->hx) * (rect->hy);
   localM[0][3] = localM[3][0] = localM[1][2] = localM[2][1] = gamma / 36. * (rect->hx) * (rect->hy); // побочная диагональ
}

void FEM::CreateG(rectangle *rect)
{
   localG[0][0] = localG[1][1] = localG[2][2] = localG[3][3] = lambda * (rect->hy / rect->hx + rect->hx / rect->hy) / 3.;
   localG[0][1] = localG[1][0] = localG[2][3] = localG[3][2] = -lambda * (rect->hy / rect->hx - rect->hx / rect->hy / 2.) / 3.;
   localG[0][2] = localG[2][0] = localG[1][3] = localG[3][1] = lambda * (rect->hy / rect->hx / 2. - rect->hx / rect->hy) / 3.;
   localG[0][3] = localG[3][0] = localG[1][2] = localG[2][1] = -lambda * (rect->hy / rect->hx + rect->hx / rect->hy) / 6.;
}

void FEM::Createb(int s) // 1/2bj + 1/2bj-2 + 2/dt^2 Mx qj-1 - 1/dt^2 Mk qj-2 + 1/2dt Mo qj-2 - 1/2 Gqj-2
{
   real localb[4];
   real f_[4];
   //for (list <rectangle*>::iterator iter = rects.begin(); iter != rects.end(); iter++)
   for (int k = 0; k < num_of_FE; k++)
   {
      rectangle* rect_ = rects[k];

      for (int i = 0; i < 4; i++)
         f_[i] = f(rect_->knots[i], s);

      localb[0] = (rect_->hx) * (rect_->hy) * (4. * f_[0] + 2. * f_[1] + 2. * f_[2] + f_[3]) / 36.;
      localb[1] = (rect_->hx) * (rect_->hy) * (2. * f_[0] + 4. * f_[1] + f_[2] + 2. * f_[3]) / 36.;
      localb[2] = (rect_->hx) * (rect_->hy) * (2. * f_[0] + f_[1] + 4. * f_[2] + 2. * f_[3]) / 36.;
      localb[3] = (rect_->hx) * (rect_->hy) * (f_[0] + 2. * f_[1] + 2. * f_[2] + 4. * f_[3]) / 36.;
      
      for (int i = 0; i < 4; i++)
         b[rect_->knots_num[i]] += localb[i];
   }
}

real FEM::scalar(real* v, real* u)
{
   real sum = 0;
   for (int i = 0; i < num_of_knots; i++)
      sum += v[i] * u[i];

   return sum;
}

void FEM::MulAb(real* v, Matrix* A, real* b)
{
   //real *out = new real[num_of_knots];

   for (int i = 0; i < num_of_knots; i++)
      v[i] = A->d[i] * b[i];

   for (int i = 0; i < num_of_knots; i++)
      for (int j = A->ig[i]; j < A->ig[i + 1]; j++) // -1?
      {
         v[i] += A->l[j] * b[A->jg[j]];
         v[A->jg[j]] += A->u[j] * b[i];
      }
}

void FEM::SolveSLAE(int s)
{
   real res, alpha, beta, skp, eps = 1e-17;
   int i, k;
   x = q[s];

   double lastres;
   MulAb(ff, A, x);
   for (i = 0; i < num_of_knots; i++)
      z[i] = r[i] = b[i] - ff[i];
   MulAb(p, A, z);
   res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));

   for (k = 1; k < 100000 && res > eps; k++)
   {
      lastres = res;
      skp = scalar(p, p);
      alpha = scalar(p, r) / skp;
      for (i = 0; i < num_of_knots; i++)
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      MulAb(ff, A, r);
      beta = -scalar(p, ff) / skp;
      for (i = 0; i < num_of_knots; i++)
      {
         z[i] = r[i] + beta * z[i];
         p[i] = ff[i] + beta * p[i];
      }
      res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));
   }
}


void FEM::Output()
{
}

real FEM::diffNorm(real* q)
{
   return real();
}

real FEM::norm(real* b)
{
   return real();
}