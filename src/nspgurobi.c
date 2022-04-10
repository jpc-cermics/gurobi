/* Nsp
 * Copyright (C) 2022-2022 Jean-Philippe Chancelier ENPC/Cermics
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

#include <gurobi_c.h>

#include <nsp/nsp.h>
#include <nsp/objects.h>
#include <nsp/imatrix.h>
#include "nspgurobi.h"

double nsp_gurobi_dbl_max(void)
{
  return  GRB_INFINITY;
}



/* creates a gurobi problem and solve the problem.
 * If a filname is given the problem is also saved 
 * if save_only if true the problem is just saved.
 */ 

int nsp_gurobi_solve(const char* problemName, int sense, int ncols, int nrows, int neq,
		    NspIMatrix*Cmatbeg, NspIMatrix *Cmatcount, NspIMatrix *Cmatind, NspMatrix *Cmatval, 
		    NspMatrix *lower, NspMatrix *upper, NspMatrix *Objective,
		    NspIMatrix *Qmatbeg,NspIMatrix *Qmatcnt, NspIMatrix *Qmatind, NspMatrix *Qmatval, 
		    NspMatrix *Rhs,const char *columnType,  NspMatrix *X,NspMatrix *Lambda,
		    NspMatrix *RetCost,NspMatrix *Retcode, char *rowType,
		    int semiCount, int *semiIndex,NspHash *Options,int loglevel,
		    const char *filename,int save_only)
{
  /* matrix A part */
  int *matrixBegin = (int *) Cmatbeg->Iv;
  int *matrixIndex = (int *) Cmatind->Iv;
  int *matrixCount = (int *) Cmatcount->Iv;
  double *matrixValues = Cmatval->R; 
  int colCount = ncols, rowCount = nrows; 
  int objectSense = (sense == 0 ) ?  GRB_MINIMIZE: GRB_MAXIMIZE;
  double *objectCoeffs = Objective->R;
  double *lowerBounds = lower->R;
  double *upperBounds = upper->R;
  double *rhsValues = Rhs->R;

  GRBenv   *env   = NULL;
  GRBmodel *model = NULL;
  int       error = 0;
  int       optimstatus;
  double    objval;
  
  double *dj=NULL,*slack=NULL;
  int ret = FAIL;

  /* Create environment */
  error = GRBloadenv(&env, NULL);
  if (error) goto QUIT;

  error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, loglevel);
  if (error) goto QUIT;

  /* Now copy the LP part of the problem data into the lp */
  error = GRBloadmodel(env, &model, "gc_pwl_c", colCount, rowCount,
                       objectSense, 0.0, objectCoeffs,
		       rowType, rhsValues,
		       matrixBegin, matrixCount, matrixIndex, matrixValues,
		       lowerBounds, upperBounds,
		       NULL, NULL, NULL);
  if (error) goto QUIT;
      
    
  /* solves the pb */
  error = GRBoptimize(model);
  if (error) goto QUIT;
  
  /* Capture solution information */
  
  error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
  if (error) goto QUIT;
  
  Retcode->R[0]= optimstatus;

  switch ( optimstatus )
    {
    case GRB_OPTIMAL:
      Sciprintf("Model is optimal\n");break;
    case GRB_INFEASIBLE:
      Sciprintf("Model is infeasible\n");break;
    case GRB_INF_OR_UNBD:
      Sciprintf("Model is infinite or unbounded\n");break;
    case GRB_UNBOUNDED:
      Sciprintf("Model is unbounded\n");break;
    }

  
  if (optimstatus == GRB_OPTIMAL)
    {
      error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
      if (error) goto QUIT;

      error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, colCount,X->R );
      if (error) goto QUIT;

      if (  columnType ) 
	{
	  int i;
	  /* no multipliers for mip */
	  for ( i = 0 ; i < rowCount; i++)
	    Lambda->R[i]=0.0;
	}
      else
	{
	  /* first the inequality constraints */
	  error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, neq, rowCount-neq, Lambda->R);
	  if ( error ) {
	    Scierror("Failed to get optimal pi values.\n");
	    goto QUIT;
	  }
	  /* then the equality constraints */
	  error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI,0, neq, Lambda->R+(rowCount-neq));
	  if ( error ) {
	    Scierror("Failed to get optimal pi values.\n");
	    goto QUIT;
	  }
	}
    
      RetCost->R[0]= objval;
    }
  else
    {
      int i;
      double d=0;d=1/d;
      RetCost->R[0]= (sense == 0 ) ? + d: - d;
      for ( i = 0 ; i < rowCount; i++) Lambda->R[i]=0.0;
      for ( i = 0 ; i < colCount; i++) X->R[i]=0.0;

    }
      
  ret = OK ;
  
 QUIT:
  /* Free up the solution */
  if ( dj != NULL) free(dj);
  if ( slack != NULL) free(slack);

  if (error) {
    Scierror("Erro: %s\n", GRBgeterrormsg(env));
    ret = FAIL;
  }

  /* Free model */
  if (model != NULL ) GRBfreemodel(model);

  /* Free environment */
  if ( env != NULL ) GRBfreeenv(env);
  
  return ret;
}
