/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   *
 *   luca.bisti@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "poly.h"
#include "tnode.h"
#include "simplex.h"

/***************************************************************
*           LINEAR PROGRAMMING: THE SIMPLEX METHOD             *
* ------------------------------------------------------------ *
* SAMPLE RUN:                                                  *
* Maximize z = x1 + x2 + 3x3 -0.5x4 with conditions:           * 
*          x1  + 2x3 <= 740                                    *
*          2x2 - 7x4 <= 0                                      *
*          x2  - x3 + 2x4 >= 0.5                               *
*          x1 + x2 + x3 +x4 = 9                                *
*          and all x's >=0.                                    *
*                                                              *
* Number of variables in E.F.: 4                               *
* Number of <= inequalities..: 2                               *
* Number of >= inequalities..: 1                               *
* Number of = equalities.....: 1                               *
* Input Economic Function:                                     *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: 3                                           *
* Coefficient # 4: -0.5                                        *
* Constant term..: 0                                           *
* Input constraint # 1:                                        *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 0                                           *
* Coefficient # 3: 2                                           *
* Coefficient # 4: 0                                           *
* Constant term..: 740                                         *
* Input constraint # 2:                                        *
* Coefficient # 1: 0                                           *
* Coefficient # 2: 2                                           *
* Coefficient # 3: 0                                           *
* Coefficient # 4: -7                                          *
* Constant term..: 0                                           *
* Input constraint # 3:                                        *
* Coefficient # 1: 0                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: -1                                          *
* Coefficient # 4: 2                                           *
* Constant term..: 0.5                                         *
* Input constraint # 4:                                        *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: 1                                           *
* Coefficient # 4: 1                                           *
* Constant term..: 9                                           *
*                                                              *
* Input Table:                                                 *
*    0.00    1.00    1.00    3.00   -0.50                      *
*  740.00   -1.00    0.00   -2.00    0.00                      *
*    0.00    0.00   -2.00    0.00    7.00                      *
*    0.50    0.00   -1.00    1.00   -2.00                      *
*    9.00   -1.00   -1.00   -1.00   -1.00                      *
*                                                              *
* Maximum of E.F. =    17.02500                                *
*  X1 =    0.000000                                            *
*  X2 =    3.325000                                            *
*  X3 =    4.725000                                            *
*  X4 =    0.950000                                            *
*                                                              *
* ------------------------------------------------------------ *
* Reference: "Numerical Recipes By W.H. Press, B. P. Flannery, *
*             S.A. Teukolsky and W.T. Vetterling, Cambridge    *
*             University Press, 1986" [BIBLI 08].              *
*                                                              *
*                       C++ Release 1.0 By J-P Moreau, Paris   *
***************************************************************/ 
#include <stdio.h>
#include <math.h>

namespace Simplex {

#define  REAL  double

typedef REAL MAT[MMAX][NMAX];

static MAT  A;
static int  IPOSV[MMAX], IZROV[NMAX];
static int  i,j,ICASE,N,M,M1,M2,M3;
static REAL R;


  void simp1(MAT,int,int *,int,int,int *,REAL *);
  void simp2(MAT,int,int,int *,int,int *,int,REAL *);
  void simp3(MAT,int,int,int,int);


  void simplx(MAT a,int m,int n,int m1,int m2,int m3,int *icase,int *izrov, int *iposv)  { 
 /*---------------------------------------------------------------------------------------- 
	  USES simp1,simp2,simp3 
	  Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3, 
	  and output parameters a, icase, izrov, and iposv are described above (see reference). 
	  Parameters: MMAX is the maximum number of constraints expected; NMAX is the maximum number 
	  of variables expected; EPS is the absolute precision, which should be adjusted to the 
	  scale of your variables.                                    
	  -----------------------------------------------------------------------------------------*/
	  int i,ip,ir,is,k,kh,kp,m12,nl1,nl2,l1[NMAX],l2[MMAX],l3[MMAX]; 
	  REAL bmax,q1,EPS=1e-6; 
	  if(m != m1+m2+m3) {
		  printf(" Bad input constraint counts in simplx.\n");
		  return;
	  }	
	  nl1=n; 
	  for (k=1; k<=n; k++) { 
		  l1[k]=k;     //Initialize index list of columns admissible for exchange.
		  izrov[k]=k;  //Initially make all variables right-hand. 
	  } 
	  nl2=m;
	  for (i=1; i<=m; i++) { 
		  if (a[i+1][1] < 0.0) {
			  printf(" Bad input tableau in simplx, Constants bi must be nonnegative.\n");
			  return;
		  }
		  l2[i]=i;
		  iposv[i]=n+i; 
/*------------------------------------------------------------------------------------------------
		  Initial left-hand variables. m1 type constraints are represented by having their slackv ariable 
		  initially left-hand, with no artificial variable. m2 type constraints have their slack 
		  variable initially left-hand, with a minus sign, and their artificial variable handled implicitly 
		  during their first exchange. m3 type constraints have their artificial variable initially 
		  left-hand.     
		  ------------------------------------------------------------------------------------------------*/ 
	  }
	  for (i=1; i<=m2; i++) l3[i]=1;
	  ir=0;
	  if(m2+m3 == 0) goto e30; //The origin is a feasible starting solution. Go to phase two. 
	  ir=1;
	  for (k=1; k<=n+1; k++) { //Compute the auxiliary objective function. 
		  q1=0.0; 
		  for (i=m1+1; i<=m; i++) q1 += a[i+1][k]; 
		  a[m+2][k]=-q1; 
	  } 
	  e10: simp1(a,m+1,l1,nl1,0,&kp,&bmax);    //Find max. coeff. of auxiliary objective fn 
	  if(bmax <= EPS && a[m+2][1] < -EPS) { 
		  *icase=-1;    //Auxiliary objective function is still negative and can’t be improved, 
		  return;       //hence no feasible solution exists.
	  }	
	  else if (bmax <= EPS && a[m+2][1] <= EPS) { 
//Auxiliary objective function is zero and can’t be improved; we have a feasible starting vector. 
//Clean out the artificial variables corresponding to any remaining equality constraints by 
//goto 1’s and then move on to phase two by goto 30.
		  m12=m1+m2+1;
		  if (m12 <= m)
			  for (ip=m12; ip<=m; ip++) 
				  if(iposv[ip] == ip+n) {   //Found an artificial variable for an equalityconstraint. 
			  simp1(a,ip,l1,nl1,1,&kp,&bmax); 
			  if(bmax > EPS) goto e1; //Exchange with column corresponding to maximum 
				  }                         //pivot element in row. 
				  ir=0;
				  m12=m12-1;
				  if (m1+1 > m12) goto e30; 
				  for (i=m1+1; i<=m1+m2; i++)     //Change sign of row for any m2 constraints 
					  if(l3[i-m1] == 1)             //still present from the initial basis. 
						  for (k=1; k<=n+1; k++) 
							  a[i+1][k] *= -1.0; 
				  goto e30;                       //Go to phase two. 
	  } 

	  simp2(a,m,n,l2,nl2,&ip,kp,&q1); //Locate a pivot element (phase one). 
                                         
	  if(ip == 0) {                         //Maximum of auxiliary objective function is 
		  *icase=-1;                          //unbounded, so no feasible solution exists.
		  return; 
	  } 
	  e1: simp3(a,m+1,n,ip,kp); 
//Exchange a left- and a right-hand variable (phase one), then update lists. 
	  if(iposv[ip] >= n+m1+m2+1) { //Exchanged out an artificial variable for an 
                               //equality constraint. Make sure it stays 
                               //out by removing it from the l1 list. 
		  for (k=1; k<=nl1; k++) 
			  if(l1[k] == kp) goto e2; 
		  e2: nl1=nl1-1; 
		  for (is=k; is<=nl1; is++)  l1[is]=l1[is+1]; 
	  }	
	  else {
		  if(iposv[ip] < n+m1+1) goto e20;
		  kh=iposv[ip]-m1-n; 
		  if(l3[kh] == 0) goto e20;  //Exchanged out an m2 type constraint. 
		  l3[kh]=0;                  //If it’s the first time, correct the pivot column 
                               //or the minus sign and the implicit 
                               //artificial variable. 
	  }   
	  a[m+2][kp+1] += 1.0; 
	  for (i=1; i<=m+2; i++)  a[i][kp+1] *= -1.0; 
	  e20: is=izrov[kp];             //Update lists of left- and right-hand variables. 
	  izrov[kp]=iposv[ip]; 
	  iposv[ip]=is; 
	  if (ir != 0) goto e10;       //if still in phase one, go back to 10. 
//End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it. 
	  e30: simp1(a,0,l1,nl1,0,&kp,&bmax); //Test the z-row for doneness. 
	  if(bmax <= EPS) {          //Done. Solution found. Return with the good news. 
		  *icase=0; 
		  return; 
	  }
	  simp2(a,m,n,l2,nl2,&ip,kp,&q1);   //Locate a pivot element (phase two). 
	  if(ip == 0) {             //Objective function is unbounded. Report and return. 
		  *icase=1; 
		  return; 
	  } 
	  simp3(a,m,n,ip,kp);       //Exchange a left- and a right-hand variable (phase two), 
	  goto e20;                 //update lists of left- and right-hand variables and 
  }                           //return for another iteration.

// The preceding routine makes use of the following utility subroutines: 

  void simp1(MAT a,int mm,int *ll,int nll,int iabf,int *kp,REAL *bmax) {
//Determines the maximum of those elements whose index is contained in the supplied list 
//ll, either with or without taking the absolute value, as flagged by iabf. 
	  int k; 
	  REAL test; 
	  *kp=ll[1]; 
	  *bmax=a[mm+1][*kp+1];
	  if (nll < 2) return; 
	  for (k=2; k<=nll; k++) { 
		  if(iabf == 0) 
			  test=a[mm+1][ll[k]+1]-(*bmax); 
		  else
			  test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax); 
		  if(test > 0.0) { 
			  *bmax=a[mm+1][ll[k]+1]; 
			  *kp=ll[k]; 
		  } 
	  } 
	  return; 
  } 

  void simp2(MAT a, int m,int n,int *l2,int nl2,int *ip,int kp,REAL *q1) { 
	  REAL EPS=1e-6; 
//Locate a pivot element, taking degeneracy into account. 
	  int i,ii,k; 
	  REAL q,q0,qp; 
	  *ip=0; 
	  if(nl2 < 1) return;
	  for (i=1; i<=nl2; i++) 
		  if (a[i+1][kp+1] < -EPS) goto e2; 
	  return;  //No possible pivots. Return with message. 
	  e2: *q1=-a[l2[i]+1][1]/a[l2[i]+1][kp+1]; 
	  *ip=l2[i];
	  if (i+1 > nl2) return; 
	  for (i=i+1; i<=nl2; i++) { 
		  ii=l2[i];
		  if(a[ii+1][kp+1] < -EPS) { 
			  q=-a[ii+1][1]/a[ii+1][kp+1]; 
			  if (q <  *q1) { 
				  *ip=ii; 
				  *q1=q; 
			  }
			  else if (q == *q1) {  //We have a degeneracy.
				  for (k=1; k<=n; k++) { 
					  qp=-a[*ip+1][k+1]/a[*ip+1][kp+1]; 
					  q0=-a[ii+1][k+1]/a[ii+1][kp+1]; 
					  if (q0 != qp) goto e6; 
				  } 
				  e6:     if (q0 < qp) *ip=ii; 
			  } 
		  } 
	  } 
	  return; 
  } 

  void simp3(MAT a,int i1,int k1,int ip,int kp) { 
//Matrix operations to exchange a left-hand and right-hand variable (see text). 
	  int ii,kk; 
	  REAL piv; 
	  piv=1.0/a[ip+1][kp+1];
	  if (i1 >= 0) 
		  for (ii=1; ii<=i1+1; ii++) 
			  if (ii-1 != ip) { 
		  a[ii][kp+1] *= piv; 
		  for (kk=1; kk<=k1+1; kk++) 
			  if (kk-1 != kp) 
				  a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1]; 
			  } 
			  for (kk=1; kk<=k1+1; kk++) 
				  if(kk-1 !=  kp) a[ip+1][kk] =-a[ip+1][kk]*piv; 
			  a[ip+1][kp+1]=piv; 
			  return; 
  } 
  
  
/***************************************************************************************************************/
  
  extern void dump_table()
  {
	  printf("Input Table:\n");
	  printf("------------------------------\n");
	  for (i=1; i<=M+1; i++) {
		  for (j=1; j<=N+1; j++)
			  printf("%8.2f", A[i][j]);
		  printf("\n");
	  }
	  printf("------------------------------\n");
  }
  
  
  /*	This is a test function which asks for problem data interactively (keyboard input)
  		It computes the MAX of the economic function
  */
  int test_simplex() {
	  printf("\n");
	  printf(" Number of variables in E.F.: "); scanf("%d", &N);
	  printf(" Number of <= inequalities..: "); scanf("%d", &M1);
	  printf(" Number of >= inequalities..: "); scanf("%d", &M2);
	  printf(" Number of = equalities.....: "); scanf("%d", &M3);

	  M=M1+M2+M3;  // Total number of constraints

	  for (i=1; i<=M+2; i++)
		  for (j=1; j<=N+1; j++)
			  A[i][j]=0.0;

	  printf(" Input Economic Function to MAXIMIZE:\n");
	  for (i=2; i<=N+1; i++) {
		  printf(" Coefficient #%d: ", i-1);
		  scanf("%lf", &A[1][i]);
	  }
	  printf(" Constant term : ");
	  scanf("%lf", &A[1][1]);
// input constraints    
	  for (i=1; i<=M; i++) {
		  printf(" Input constraint #%d: \n", i);
		  for (j=2; j<=N+1; j++) {
			  printf(" Coefficient #%d: ", j-1); 
			  scanf("%lf", &R);
			  A[i+1][j] = -R;
		  }
		  printf(" Constant term : ");
		  scanf("%lf", &A[i+1][1]);   
	  }

	  printf("\n Input Table:\n");
	  for (i=1; i<=M+1; i++) {
		  for (j=1; j<=N+1; j++)
			  printf("%8.2f", A[i][j]);
		  printf("\n");
	  }

	  simplx(A,M,N,M1,M2,M3,&ICASE,IZROV,IPOSV);

	  if (ICASE==0) {  //result ok.
		  printf("\n Maximum of E.F. = %f\n", A[1][1]);

		  for (i=1; i<=N; i++) {
			  for (j=1; j<=M; j++)
				  if (IPOSV[j] == i)  {
				  printf("  X%d = %f\n", i, A[j+1][1]);
				  goto e3;
				  }
				  printf("  X%d = %f\n", i, 0.0);
				  e3:;}
	  }
	  else
		  printf(" No solution (error code = %d).\n", ICASE);
  
	  printf("\n");
	  return 0;
  }
  
  
  /*	This method runs the simplex algorithm to minimize the passed economic function
  		EcFun = Economic Function as returned by the TNode class (the first variable will be skipped)
  		c_set = set of constraints in the form ax+b<=0 or ax+b<=0     (x>=0 are assumed implicitly)
  		dead_var = the variable at this index must be excluded from the simplex
  		MinValue = the optimum value of the econ.function will be stored here
  		Solution = array of doubles were the optimal solution X will be stored (can be NULL)
  */
  int minimum(Poly &EcFun, ConstraintSet &c_set, const uint dead_var, Poly &Optimum)
  {
	  int result = 0;
	  int nvars = EcFun.getNumVariables();
	  uint idx;
	  int num_constraints;
	  ConstraintSet gte;
	  ConstraintSet lte;
	  double solution[nvars];
	  
	  if(nvars > NMAX)
	  {
		  printf("Simplex error: max number of variables exceeded (%i, max %i)\n",nvars,NMAX);
		  return SIMPLEX_RESULT_UNFEASIBLE;
	  }
	  
	  //Optimum.Zero();
	  num_constraints = c_set.size();
	  if(num_constraints > MMAX)
	  {
		  printf("Simplex error: max number of constraints exceeded (N=%i)\n", num_constraints);
		  return SIMPLEX_RESULT_UNFEASIBLE;
	  }
	  
	  // split the constraint set into two distinct sets (for <= and >= types)
	  for(i=0; i<c_set.size(); i++)
	  {
		  CPoly c1(c_set[i]);
		  if(-c1.GetB() < 0.0) c1.InvertDirection(); // negate B because the library expects ax<=b or ax>=b
		  if(c1.greater_than_zero) gte.push_back(c1);
		  else lte.push_back(c1);
	  }
	  
	  N = nvars - 1; // discard the dead variable
	  M1 = lte.size(); // number of <= constraints
	  M2 = gte.size(); // number of >= constraints
	  M3 = 0; // number of == constraints
	  M = M1 + M2 + M3;  // total number of constraints
	  
	  // reset problem data
	  for (i=1; i<=M+2; i++)
		  for (j=1; j<=N+1; j++)
			  A[i][j]=0.0;
	  
	  // store economic function data (negate all coefficients, since min(cx) = -max(-cx)
	  A[1][1] = -EcFun.GetB(); // constant term
	  idx=2;
	  for (i=0; i<N+1; i++)
	  {
		  if(i == dead_var) continue;
		  A[1][idx] = -EcFun.GetCoeff(i);
		  idx++;
	  }
	  
	  // store constraints data
	  for (i=1; i<=M; i++)
	  {
		  CPoly &cp = (i<=M1) ? lte[i-1] : gte[i-M1-1];
		  //ith constraint
		  A[i+1][1] = -cp.GetB(); // constant term: negate it because the library expects ax<=b or ax>=b
		  idx=2;
		  for (j=0; j<N+1; j++)
		  {
			  if(j == dead_var) continue;
			  R = cp.GetCoeff(j);
			  A[i+1][idx] = -R;
			  idx++;
		  }
	  }
	  
	  //dump_table();
	  
	  // run the simplex!
	  simplx(A,M,N,M1,M2,M3,&ICASE,IZROV,IPOSV);
	  
	  if(ICASE == 0)
	  {
		  // simplex was feasible
		  Optimum.SetB(-A[1][1]);

		  for (i=1; i<=N; i++) {
			  for (j=1; j<=M; j++)
				  if (IPOSV[j] == i)  {
				  solution[i-1] = A[j+1][1];
				  goto e3;
				  }
				  solution[i-1] = 0.0;
				  e3:;}

		  j=0;
		  for(i=0; i<Optimum.getNumVariables(); i++)
		  {
			  if(i == dead_var) continue;
			  Optimum.SetCoeff(i, solution[j]);
			  j++;
		  }
		  //printf("SOLUTION = "); Optimum.Print(); printf("\n");
		  result = SIMPLEX_RESULT_OK;
	  }
	  else {
		  // simplex was unfeasible
		  for(i=0; i<N; i++) solution[i]=0.0;
		  if(ICASE == 1) result = SIMPLEX_RESULT_UNBOUNDED;
		  else if (ICASE == -1) result = SIMPLEX_RESULT_UNFEASIBLE;
		  else result = SIMPLEX_RESULT_UNK_ERROR;
	  }
	  
	  return result;
  }
  
  
  void test_es4()
  {
	  uint nvars = 5;
	  int dead_var = 0;
	  Poly fobj(nvars);
	  CPoly c1(nvars);
	  Poly opt(nvars);
	  ConstraintSet cset;
	  bool result;
	  double minvalue = 0.0;
	  
	  fobj.SetCoeff(0, 0.0); //dead variable!
	  fobj.SetCoeff(1, 0.1);
	  fobj.SetCoeff(2, 0.9);
	  fobj.SetCoeff(3, 0.08);
	  fobj.SetCoeff(4, 1.2);
	  fobj.SetB(0.0);
	  
	  c1.SetCoeff(1, 9.9);
	  c1.SetCoeff(2, 20.0);
	  c1.SetCoeff(3, 2.0);
	  c1.SetCoeff(4, 22.0);
	  c1.SetB(70.0);
	  cset.push_back(c1);
	  
	  c1.SetCoeff(1, 1.2);
	  c1.SetCoeff(2, 1.1);
	  c1.SetCoeff(3, 0.2);
	  c1.SetCoeff(4, 26.0);
	  c1.SetB(50.0);
	  cset.push_back(c1);
	  
	  c1.SetCoeff(1, 75.3);
	  c1.SetCoeff(2, 0.0);
	  c1.SetCoeff(3, 4.0);
	  c1.SetCoeff(4, 0.0);
	  c1.SetB(250.0);
	  cset.push_back(c1);
	  
	  c1.SetCoeff(1, 0.0);
	  c1.SetCoeff(2, 0.0);
	  c1.SetCoeff(3, 1.0);
	  c1.SetCoeff(4, 0.0);
	  c1.SetB(5.0);
	  c1.greater_than_zero = false;
	  cset.push_back(c1);
	  
	  printf("*** TEST: ES4 PROBLEM\nSimplex ");
	  result = minimum(fobj, cset, dead_var, opt);
	  
	  if(result)
	  {
		  printf("has computed the optimum: FOBJ = %lf\n", opt.GetB());
	  }
	  else printf("did not find a solution\n");
	  
  }
  
} //end of namespace Simplex
