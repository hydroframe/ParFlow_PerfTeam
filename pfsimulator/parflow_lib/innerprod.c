/*BHEADER*********************************************************************
 *
 *  Copyright (c) 1995-2009, Lawrence Livermore National Security,
 *  LLC. Produced at the Lawrence Livermore National Laboratory. Written
 *  by the Parflow Team (see the CONTRIBUTORS file)
 *  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
 *
 *  This file is part of Parflow. For details, see
 *  http://www.llnl.gov/casc/parflow
 *
 *  Please read the COPYRIGHT file or Our Notice and the LICENSE file
 *  for the GNU Lesser General Public License.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (as published
 *  by the Free Software Foundation) version 2.1 dated February 1999.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
 *  and conditions of the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA
 **********************************************************************EHEADER*/
/*****************************************************************************
* Inner Product of two vectors
*
*****************************************************************************/

#include "parflow_config.h"

#ifdef USING_PARALLEL
extern "C"{
#endif

#include "parflow.h"
#include "pf_parallel.h"


double   InnerProd(
                   Vector *x,
                   Vector *y)
{
  Grid         *grid = VectorGrid(x);
  Subgrid      *subgrid;

  Subvector    *y_sub;
  Subvector    *x_sub;

  double       *yp, *xp;

  /*
    @MCB:
    For OpenMP parallelism to work correctly here when
    calling InnerProd from an outer parallel region,
    result needs to be declared as a static variable and
    reset each time.  This allows us to use a reduction
    clause.

    WARNING: Because this is static, any kind of parallel tasking
    calls CANNOT be done if they make a call to InnerProd. Any call
    to InnerProd MUST be done by one active region at a time.
  */
  static double result = 0.0;
  double local_result = 0.0;

  int ix, iy, iz;
  int nx, ny, nz;
  int nx_v, ny_v, nz_v;

  int i_s, i, j, k, iv;

  amps_Invoice result_invoice;


  MASTER(result_invoice = amps_NewInvoice("%d", &result));

  ForSubgridI(i_s, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, i_s);

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);
    iz = SubgridIZ(subgrid);

    nx = SubgridNX(subgrid);
    ny = SubgridNY(subgrid);
    nz = SubgridNZ(subgrid);

    y_sub = VectorSubvector(y, i_s);
    x_sub = VectorSubvector(x, i_s);

    nx_v = SubvectorNX(y_sub);
    ny_v = SubvectorNY(y_sub);
    nz_v = SubvectorNZ(y_sub);

    yp = SubvectorElt(y_sub, ix, iy, iz);
    xp = SubvectorElt(x_sub, ix, iy, iz);

    iv = 0;
    // @MCB: Change this to atomic writes, its (supposedly) faster and more accurate
    #if 0
    __BoxLoopReduceI1(NO_LOCALS, result,
                      i, j, k, ix, iy, iz, nx, ny, nz,
                      iv, nx_v, ny_v, nz_v, 1, 1, 1,
    {
      result += yp[iv] * xp[iv];
    });
    #else
    _BoxLoopI1(NoWait, NO_LOCALS,
               i, j, k, ix, iy, iz, nx, ny, nz,
               iv, nx_v, ny_v, nz_v, 1, 1, 1,
    {
        local_result += yp[iv] * xp[iv];
    });

    #pragma omp atomic update
    result += local_result;

    #endif
  }

  BARRIER;

  MASTER(
  amps_AllReduce(amps_CommWorld, result_invoice, amps_Add);
  amps_FreeInvoice(result_invoice);
  IncFLOPCount(2 * VectorSize(x) - 1);
    );

  BARRIER;

  /*
    @MCB:
    This is basically to confuse the compiler so that
    it doesn't optimize away resetting result to 0
    for subsequent calls
  */
  double temp = result;
  result = 0.0;
  return temp;
}

#ifdef USING_PARALLEL
} // Extern C
#endif
