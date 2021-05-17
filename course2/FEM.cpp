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

   rectangle* rect;
   frects >> num_of_FE;
   rects = new rectangle * [num_of_FE];
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
      fbounds >> num;	//íîìåð óñëîâèÿ
      cond->bound_num = num;
      int number;
      for (int i = 0; i < 2; i++)
      {
         fbounds >> number;									// 
         cond->knot_nums[i] = number - 1;				// íîìåðà óçëîâ ðåáðà, íà êîòîðîì êðàåâîå...
         cond->knots[i] = knots[number - 1];	// 
      }
      if (cond->bound_num == 3)
         fbounds >> cond->bound_param;									//
      bounds.push_back(cond);
   }
   fbounds.close();

   fparams >> sigma >> khi >> un >> nt;
   t = new real[nt];
   for (int i = 0; i < nt; i++)
      fparams >> t[i];
   fparams.close();
#pragma endregion

   MakeSparseFormat();
   q = new real[num_of_knots]{};
   q0 = new real[num_of_knots]{};
   q1 = new real[num_of_knots]{};
   b = new real[num_of_knots]{};
   b0 = new real[num_of_knots]{};
   temp = new real[num_of_knots]{};
   ff = new real[num_of_knots]{};
   z = new real[num_of_knots]{};
   r = new real[num_of_knots]{};
   p = new real[num_of_knots]{};
}

void FEM::MakeSparseFormat()
{
   const int N = 4;
   int* list1, * list2;
   int* listbeg = new int[num_of_knots];

   for (int i = 0; i < num_of_knots; i++)
      listbeg[i] = -1;

   list1 = new int[num_of_knots * num_of_knots]{};
   list2 = new int[num_of_knots * num_of_knots]{};
   int iaddr, ind1, ind2, k;
   num_of_not_zero = -1;

   for (int iel = 0; iel < num_of_FE; iel++) // ïðîõîä ïî âñåì êý
   {
      for (int i = 0; i < 4; i++) // ïî óçëàì êý
      {
         k = rects[iel]->knots_num[i]; //
         for (int j = i + 1; j < 4; j++) // 
         {
            ind1 = k;
            ind2 = rects[iel]->knots_num[j];  //
            if (ind2 < ind1) //çàíåñåíèå ñâÿçè áîëüøîãî íîìåðà ñ ìåíüøèì, òå ind2 ñ ind1
            {
               ind1 = ind2;
               ind2 = k;
            }
            iaddr = listbeg[ind2];
            if (iaddr == -1) // åñëè ñïèñîê ïóñò
            {
               num_of_not_zero++;
               listbeg[ind2] = num_of_not_zero;
               list1[num_of_not_zero] = ind1;
               list2[num_of_not_zero] = -1;
            }
            else // èùåì â ñïèñêå ind1
            {
               while (list1[iaddr] < ind1 && list2[iaddr] >= 0)
                  iaddr = list2[iaddr];
               if (list1[iaddr] > ind1)  // åñëè íå íàøëè è âñòðåòèëè ýë ñ áîëüøèì íîìåðîì 
               {                         // äîáàâëÿåì ïåðåä íèì, ÷òîá ñïèñîê áûë óïîðÿäî÷åííûì
                  num_of_not_zero++;
                  list1[num_of_not_zero] = list1[iaddr];
                  list2[num_of_not_zero] = list2[iaddr];
                  list1[iaddr] = ind1;
                  list2[iaddr] = num_of_not_zero;
               }
               else if (list1[iaddr] < ind1) // åñëè íå íàøëè, òî ïèõàåì â êîíåö
               {
                  num_of_not_zero++;
                  list2[iaddr] = num_of_not_zero;
                  list1[num_of_not_zero] = ind1;
                  list2[num_of_not_zero] = -1;
               }
            }
         }
      }
   }
   num_of_not_zero++;
   A = new Matrix;
   G = new Matrix;
   M = new Matrix;
   G->ig = M->ig = A->ig = new int[num_of_knots + 1]{};
   G->jg = M->jg = A->jg = new int[num_of_not_zero]{};  // +1???

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
   for (int i = 0; i < num_of_not_zero; i++)
      cout << A->jg[i] << " ";
   cout << "\nig: ";
   for (int i = 0; i < num_of_knots + 1; i++)
      cout << A->ig[i] << " ";
   cout << "\n";

   A->l = new real[num_of_not_zero]{};
   M->l = new real[num_of_not_zero]{};
   G->l = new real[num_of_not_zero]{};

   A->u = new real[num_of_not_zero]{};
   M->u = new real[num_of_not_zero]{};
   G->u = new real[num_of_not_zero]{};

   A->d = new real[num_of_knots]{};
   M->d = new real[num_of_knots]{};
   G->d = new real[num_of_knots]{};

}

void FEM::AddElement(Matrix* A, int i, int j, real elem)
{
   int k = 0;
   if (i == j)
      A->d[i] += elem;
   else if (i < j)
   {
      for (k = A->ig[j]; k < A->ig[j + 1]; k++)
         if (A->jg[k] == i) break;

      A->u[k] += elem; // i-1?
   }
   else
   {
      for (k = A->ig[i]; k < A->ig[i + 1]; k++)
         if (A->jg[k] == j) break;

      A->l[k] += elem; // i-1??
   }

}

void FEM::AddLocal(Matrix* A, int knot_num[4], real localA[4][4], real coeff)
{
   int ibeg, iend, ind;
   for (int i = 0; i < 4; i++)
      A->d[knot_num[i]] += coeff * localA[i][i];
   for (int i = 0; i < 4; i++)
   {
      ibeg = A->ig[knot_num[i]];
      for (int j = 0; j < i; j++) // i - 1?
      {
         iend = A->ig[knot_num[i] + 1] - 1;  // -1 ?
         while (A->jg[ibeg] != knot_num[j])
         {
            ind = (ibeg + iend) / 2 + 1;
            if (A->jg[ind] <= knot_num[j])
               ibeg = ind;
            else
               iend = (iend == ind) ? ind - 1 : ind;
         }
         A->l[ibeg] += coeff * localA[i][j];
         A->u[ibeg] += coeff * localA[j][i];
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
   ofstream out("Result.txt");
   out.scientific;
   out.precision(14);
   cout.scientific;
   cout.precision(14);

   char title[] = "\n| q*\t\t\t| q\t\t\t| |q*-q|\t\t|\n+-----------------------+-----------------------+-----------------------+\n";
   out << title;
   cout << title;
#pragma region initial_q
   for (int i = 0; i < num_of_knots; i++)
      q[i] = q0[i] = ug(knots[i], t[0]);
   Output(0, out);
   for (int i = 0; i < num_of_knots; i++)
      q[i] = q1[i] = ug(knots[i], t[1]);
   Output(1, out);
#pragma endregion

   for (int s = 2; s < nt; s++)
   {
      CreateA(s);
      Createb(s);
      MakeBounds(s);
      SolveSLAE(s);
      Output(s, out);

      equalize(q0, q1);
      equalize(q1, q);
   }

   out.close();
}

void FEM::Output(int s, ofstream& out)
{
   cout << "t = " << t[s] << '\n';
   out << "t = " << t[s] << '\n';
   for (int i = 0; i < num_of_knots; i++)
   {
      out << scientific << "| " << ug(knots[i], t[s]) << "\t| " << q[i] << "\t| "
         << abs(q[i] - ug(knots[i], t[s])) << "\t|\n";
      cout << scientific << "| " << ug(knots[i], t[s]) << "\t| " << q[i] << "\t| "
         << abs(q[i] - ug(knots[i], t[s])) << "\t|\n";
   }
   cout << "+-----------------------+-----------------------+-----------------------+\n";
   out << "+-----------------------+-----------------------+-----------------------+\n";

}

void FEM::MakeBounds(int s)
{
   for (list<bound*>::iterator iter = bounds.begin(); iter != bounds.end(); iter++)
   {
      bound* cond = *iter;

      if (cond->bound_num == 1) // 1 краевое
      {
         for (int i = 0; i < 2; i++)
         {
            A->d[cond->knot_nums[i]] = 1;
            for (int j = A->ig[cond->knot_nums[i]]; j < A->ig[cond->knot_nums[i] + 1]; j++)
               A->l[j] = 0.;
            for (int j = 0; j < A->ig[num_of_knots]; j++)
               if (A->jg[j] == cond->knot_nums[i])
                  A->u[j] = 0.;

            b[cond->knot_nums[i]] = ug(cond->knots[i], t[s]);
         }
      }
      else if(cond->bound_num == 2)
      {
         int h = abs(cond->knot_nums[0] - cond->knot_nums[0]) != 1 ? abs(cond->knots[0].x - cond->knots[1].x) : abs(cond->knots[0].y - cond->knots[1].y);
         real b_add[2] = { h / 6. * (2. * theta(cond->knots[0], t[s]) + theta(cond->knots[1], t[s])),
                           h / 6. * (theta(cond->knots[0], t[s]) + 2. * theta(cond->knots[1], t[s])) };
         b[cond->knot_nums[0]] += b_add[0];
         b[cond->knot_nums[1]] += b_add[1];
      }

      else if (cond->bound_num == 3)
      {
         int h = cond->bound_param * ( abs(cond->knot_nums[0] - cond->knot_nums[0]) != 1 ? abs(cond->knots[0].x - cond->knots[1].x) : abs(cond->knots[0].y - cond->knots[1].y)) / 6.;
         real b_add[2] = {  h * (2 * ubeta(cond->knots[0], t[s]) + ubeta(cond->knots[1], t[s])),
                            h * (ubeta(cond->knots[0], t[s]) + 2 * ubeta(cond->knots[1], t[s])) };
         b[cond->knot_nums[0]] += b_add[0];
         b[cond->knot_nums[1]] += b_add[1];

         real localAS3[2][2] = {{ 2. * h , h}, { h, 2. * h }};

         for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
               AddElement(A, i, j, localAS3[i][j]);
      }
   }
}

void FEM::CreateA(int s) // 1/dt^2 Mx + 1/2dt Mo + 1/2 G
{
   rectangle* rect;
   real dt = (t[s] - t[s - 1]);

   for (int i = 0; i < num_of_knots; i++)
      G->d[i] = A->d[i] = 0.;
   for (int i = 0; i < num_of_not_zero; i++)
      G->u[i] = G->l[i] = A->l[i] = A->u[i] = 0.;

   for (int i = 0; i < num_of_FE; i++)
   {
      rect = rects[i];
      CreateG(rect);
      CreateM(rect);
      AddLocal(A, rect->knots_num, localM, khi / (dt * dt) + sigma / (2. * dt)); // A = 1/dt^2 Mx + 1/ 2dt Mo
      AddLocal(G, rect->knots_num, localG, .5);
   }
   for (int i = 0; i < num_of_not_zero; i++) // A += 1/2 G
   {
      A->l[i] += G->l[i];
      A->u[i] += G->u[i];
   }
   for (int i = 0; i < num_of_knots; i++)
      A->d[i] += G->d[i];
}

void FEM::CreateM(rectangle* rect)
{
   localM[0][0] = localM[1][1] = localM[2][2] = localM[3][3] = (rect->hx) * (rect->hy) / 9.;
   localM[0][1] = localM[1][0] = localM[0][2] = localM[2][0] =
   localM[2][3] = localM[3][2] = localM[1][3] = localM[3][1] = (rect->hx) * (rect->hy) / 18.;
   localM[0][3] = localM[3][0] = localM[1][2] = localM[2][1] = (rect->hx) * (rect->hy) / 36.;
}

void FEM::CreateG(rectangle* rect) // надо сделать разложение лямбды
{
   localG[0][0] = localG[1][1] = localG[2][2] = localG[3][3] = lambda(rect->knots[0], 0) * (rect->hy / rect->hx + rect->hx / rect->hy) / 3.;
   localG[0][1] = localG[1][0] = localG[2][3] = localG[3][2] = -lambda(rect->knots[1], 0) * (rect->hy / rect->hx - rect->hx / rect->hy / 2.) / 3.;
   localG[0][2] = localG[2][0] = localG[1][3] = localG[3][1] = lambda(rect->knots[2], 0) * (rect->hy / rect->hx / 2. - rect->hx / rect->hy) / 3.;
   localG[0][3] = localG[3][0] = localG[1][2] = localG[2][1] = -lambda(rect->knots[3], 0) * (rect->hy / rect->hx + rect->hx / rect->hy) / 6.;
}

void FEM::Createb(int s) // 1/2bj + 1/2bj-2 + 2/dt^2 Mx qj-1 - 1/dt^2 Mk qj-2 + 1/2dt Mo qj-2 - 1/2 Gqj-2
{
   real localb[4];
   real d[4];    // M ((bj + bj-2)/2 + 2x/dt^2 qj-1 - x/dt^2 qj-2 + o/2dt qj-2) // тк o и x - константы
   real dt = (t[s] - t[s - 1]);
   nullify(temp);
   nullify(b);
   MulAb(temp, G, q0);  // 1/2 G qj-2

   for (int k = 0; k < num_of_FE; k++)
   {
      rectangle* rect_ = rects[k];

      for (int i = 0; i < 4; i++)
         d[i] = (f(rect_->knots[i], t[s]) + f(rect_->knots[i], t[s - 2])) / 2.                  // (fj + fj+2) / 2 
              + (khi / (dt * dt)) * (2. * q1[rect_->knots_num[i]] - q0[rect_->knots_num[i]])    // 2/dt^2 Mx qj-1 - 1/dt^2 Mk qj-2
              + sigma * q0[rect_->knots_num[i]] / (2. * dt);                                    // 1/2dt Mo qj-2
      
      localb[0] = ((rect_->hx) * (rect_->hy) * (4. * d[0] + 2. * d[1] + 2. * d[2] + d[3]) / 36.) - temp[rect_->knots_num[0]] ;
      localb[1] = ((rect_->hx) * (rect_->hy) * (2. * d[0] + 4. * d[1] + d[2] + 2. * d[3]) / 36.) - temp[rect_->knots_num[1]] ;
      localb[2] = ((rect_->hx) * (rect_->hy) * (2. * d[0] + d[1] + 4. * d[2] + 2. * d[3]) / 36.) - temp[rect_->knots_num[2]] ;
      localb[3] = ((rect_->hx) * (rect_->hy) * (d[0] + 2. * d[1] + 2. * d[2] + 4. * d[3]) / 36.) - temp[rect_->knots_num[3]] ;

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

void FEM::nullify(real* x)
{
   for (int i = 0; i < num_of_knots; i++)
      x[i] = 0.;
}

void FEM::SolveSLAE(int s)
{
   real res, alpha, beta, skp, eps = 1e-167;
   int i, k;

   MulAb(ff, A, q);
   for (i = 0; i < num_of_knots; i++)
      z[i] = r[i] = b[i] - ff[i];
   MulAb(p, A, z);
   res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));

   for (k = 1; k < 100000 && res > eps; k++)
   {
      skp = scalar(p, p);
      alpha = scalar(p, r) / skp;
      for (i = 0; i < num_of_knots; i++)
      {
         q[i] += alpha * z[i];
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

      if (k % 100)
      {
          MulAb(ff, A, q);
          for (i = 0; i < num_of_knots; i++)
              z[i] = r[i] = b[i] - ff[i];
          MulAb(p, A, z);
          res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));
      }
   }
}
