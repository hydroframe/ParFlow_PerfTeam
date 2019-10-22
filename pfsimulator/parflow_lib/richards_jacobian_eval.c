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

/*
 * SGS TODO this needs some work in the overland flow current
 * implemnetation is doing communication and computations that are not
 * needed.  The C matrix is wasting a lot of space and communication.
 * Making it a set of vectors will save a great deal of space or adding
 * 2D Matrix support.
 *
 * SGS TODO There is a problem attempting to avoid doing the overland
 * flow additions since the flag is not local.  A neighbor doing
 * overland flow means the process does as well if the overland flow
 * cell is on the boundary.
 */

#include "parflow.h"
#include "llnlmath.h"
#include "llnltyps.h"
#include "assert.h"

#include "timer.h"
#include "get_modules_xtra.h"
#include "inline_defs.h"
/*---------------------------------------------------------------------
 * Define module structures
 *---------------------------------------------------------------------*/

// Which Jacobian to use.
//
enum JacobianType {
  no_nonlinear_jacobian,
  not_set,
  simple,
  overland_flow
};

typedef struct {
  enum JacobianType type;
  double SpinupDampP1; // NBE
  double SpinupDampP2; // NBE
  int tfgupwind;  // @RMM
} PublicXtra;

typedef struct {
  Problem      *problem;

  PFModule     *density_module;
  PFModule     *saturation_module;
  PFModule     *rel_perm_module;
  PFModule     *bc_pressure;
  PFModule     *bc_internal;
  PFModule     *overlandflow_module;  //DOK
  PFModule     *overlandflow_module_diff;  //@LEC
  PFModule     *overlandflow_module_kin;

  /* The analytic Jacobian matrix is decomposed as follows:
   *
   *      [ JC  JE ]
   * J =  |        |
   *      [ JF  JB ]
   *
   * where JC corresponds to surface-surface interations,
   *      JB corresponds to subsurface-subsurface interations,
   *      JE corresponds to surface-subsurface interations, and
   *      JF corresponds to subsurface-surface interations.
   *
   * To make for a more efficient implementation, we store the
   * interactions for JE and JF as part of JB, so that JC handles
   * only surface-surface interactions, and JB handles the rest
   *
   * To make this more general, JB = J whenever there is no
   * overland flow contributions to the Jacobian. Hence the
   * analytic Jacobian for the subsurface flow is invoked instead.
   */
  Matrix       *J;
  Matrix       *JC;

  Grid         *grid;
  double       *temp_data;
} InstanceXtra;

/*--------------------------------------------------------------------------
 * Static stencil-shape definition
 *--------------------------------------------------------------------------*/

int jacobian_stencil_shape[7][3] = { { 0, 0, 0 },
                                     { -1, 0, 0 },
                                     { 1, 0, 0 },
                                     { 0, -1, 0 },
                                     { 0, 1, 0 },
                                     { 0, 0, -1 },
                                     { 0, 0, 1 } };


int jacobian_stencil_shape_C[5][3] = { { 0, 0, 0 },
                                       { -1, 0, 0 },
                                       { 1, 0, 0 },
                                       { 0, -1, 0 },
                                       { 0, 1, 0 } };

/*---------------------------------------------------------------------
 * Define macros for jacobian evaluation
 *---------------------------------------------------------------------*/

#define PMean(a, b, c, d)    HarmonicMean(c, d)
#define PMeanDZ(a, b, c, d)     HarmonicMeanDZ(a, b, c, d)
#define RPMean(a, b, c, d)   UpstreamMean(a, b, c, d)
#define Mean(a, b) ArithmeticMean(a, b)  //@RMM

/*  This routine provides the interface between KINSOL and ParFlow
 *  for richards' equation jacobian evaluations and matrix-vector multiplies.*/

int       KINSolMatVec(
                       void *   current_state,
                       N_Vector x,
                       N_Vector y,
                       int *    recompute,
                       N_Vector pressure)
{
  PFModule    *richards_jacobian_eval = StateJacEval(((State*)current_state));
  Matrix      *J = StateJac(((State*)current_state));
  Matrix      *JC = StateJacC(((State*)current_state));
  Vector      *old_pressure = StateOldPressure(((State*)current_state));
  Vector      *saturation = StateSaturation(((State*)current_state));
  Vector      *density = StateDensity(((State*)current_state));
  ProblemData *problem_data = StateProblemData(((State*)current_state));
  double dt = StateDt(((State*)current_state));
  double time = StateTime(((State*)current_state));

  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(richards_jacobian_eval);

  PFModule    *bc_pressure = (instance_xtra->bc_pressure);

  StateBCPressure((State*)current_state) = bc_pressure;

  InitVector(y, 0.0);

  /*
   * Compute Jacobian if needed.
   */
  if (*recompute)
  {
    /*
      @MCB:
      This should be a direct function call, should it not?
      There's only one function with this typedef.
    */
    PFModuleInvokeType(RichardsJacobianEvalInvoke, richards_jacobian_eval,
                       (pressure, old_pressure, &J, &JC, saturation, density, problem_data,
                        dt, time, 0));

    *recompute = 0;
    StateJac(((State*)current_state)) = J;
    StateJacC(((State*)current_state)) = JC;
  }

  if (JC == NULL)
    Matvec(1.0, J, x, 0.0, y);
  else
    MatvecSubMat(current_state, 1.0, J, JC, x, 0.0, y);

  return(0);
}


/*  This routine evaluates the Richards jacobian based on the current
 *  pressure values.  */

void    RichardsJacobianEval(
                             Vector *     pressure, /* Current pressure values */
                             Vector *     old_pressure, /* Pressure values at previous timestep */
                             Matrix **    ptr_to_J, /* Pointer to the J pointer - this will be set
                                                     * to instance_xtra pointer at end */
                             Matrix **    ptr_to_JC, /* Pointer to the JC pointer - this will be set
                                                      * to instance_xtra pointer at end */
                             Vector *     saturation, /* Saturation / work vector */
                             Vector *     density, /* Density vector */
                             ProblemData *problem_data, /* Geometry data for problem */
                             double       dt, /* Time step size */
                             double       time, /* New time value */
                             int          symm_part) /* Specifies whether to compute just the
                                                      * symmetric part of the Jacobian (1), or the
                                                      * full Jacobian */
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  Problem     *problem = (instance_xtra->problem);

  PFModule    *density_module = (instance_xtra->density_module);
  PFModule    *saturation_module = (instance_xtra->saturation_module);
  PFModule    *rel_perm_module = (instance_xtra->rel_perm_module);
  PFModule    *bc_pressure = (instance_xtra->bc_pressure);
  PFModule    *bc_internal = (instance_xtra->bc_internal);
  PFModule    *overlandflow_module = (instance_xtra->overlandflow_module);
  PFModule    *overlandflow_module_diff = (instance_xtra->overlandflow_module_diff);
  PFModule    *overlandflow_module_kin = (instance_xtra->overlandflow_module_kin);

  Matrix      *J = (instance_xtra->J);
  Matrix      *JC = (instance_xtra->JC);

  Vector      *density_der = NULL;
  Vector      *saturation_der = NULL;

  /* Re-use vectors to save memory */
  Vector      *rel_perm = NULL;
  Vector      *rel_perm_der = NULL;

  Vector      *porosity = ProblemDataPorosity(problem_data);
  Vector      *permeability_x = ProblemDataPermeabilityX(problem_data);
  Vector      *permeability_y = ProblemDataPermeabilityY(problem_data);
  Vector      *permeability_z = ProblemDataPermeabilityZ(problem_data);
  Vector      *sstorage = ProblemDataSpecificStorage(problem_data);           //sk
  Vector      *top = ProblemDataIndexOfDomainTop(problem_data);               //DOK
  Vector      *slope_x = ProblemDataTSlopeX(problem_data);                //DOK

  /* Overland flow variables */  //DOK
  Vector      *KW, *KE, *KN, *KS, *KWns, *KEns, *KNns, *KSns;
  Subvector   *kw_sub, *ke_sub, *kn_sub, *ks_sub, *kwns_sub, *kens_sub, *knns_sub, *ksns_sub, *top_sub, *sx_sub;
  double      *kw_der, *ke_der, *kn_der, *ks_der, *kwns_der, *kens_der, *knns_der, *ksns_der;

  double gravity = ProblemGravity(problem);
  double viscosity = ProblemPhaseViscosity(problem, 0);

  /* @RMM terrain following grid slope variables */
  Vector      *x_ssl = ProblemDataSSlopeX(problem_data);               //@RMM
  Vector      *y_ssl = ProblemDataSSlopeY(problem_data);               //@RMM
  Subvector   *x_ssl_sub, *y_ssl_sub;    //@RMM
  double      *x_ssl_dat, *y_ssl_dat;     //@RMM

  /* @RMM variable dz multiplier */
  Vector      *z_mult = ProblemDataZmult(problem_data);              //@RMM
  Subvector   *z_mult_sub;    //@RMM
  double      *z_mult_dat;    //@RMM

  /* @RMM Flow Barrier / Boundary values */
  Vector      *FBx = ProblemDataFBx(problem_data);
  Vector      *FBy = ProblemDataFBy(problem_data);
  Vector      *FBz = ProblemDataFBz(problem_data);
  Subvector   *FBx_sub, *FBy_sub, *FBz_sub;    //@RMM
  double      *FBx_dat, *FBy_dat, *FBz_dat;     //@RMM

  Subgrid     *subgrid;

  Subvector   *p_sub, *d_sub, *s_sub, *po_sub, *rp_sub, *ss_sub;
  Subvector   *permx_sub, *permy_sub, *permz_sub, *dd_sub, *sd_sub, *rpd_sub;
  Submatrix   *J_sub;
  Submatrix   *JC_sub;

  Grid        *grid = VectorGrid(pressure);
  Grid        *grid2d = VectorGrid(slope_x);

  double      *pp, *sp, *sdp, *pop, *dp, *ddp, *rpp, *rpdp;
  double      *permxp, *permyp, *permzp;
  double      *cp, *wp, *ep, *sop, *np, *lp, *up, *op = NULL, *ss;

  double      *cp_c, *wp_c, *ep_c, *sop_c, *np_c, *top_dat;  //DOK

  int i, j, k, r, is;
  int ix, iy, iz;
  int nx, ny, nz;
  int nx_v, ny_v, nz_v;
  int nx_m, ny_m, nz_m;
  int nx_po, ny_po, nz_po;
  int sy_v, sz_v;
  int sy_m, sz_m;
  int ip, ipo, im, iv;

  // Phase density iterators
  int id;
  int nx_d, ny_d, nz_d;
  int nx_p, ny_p, nz_p;

  int diffusive;             //@LEC

  diffusive = GetIntDefault("OverlandFlowDiffusive", 0);

  int overlandspinup;              //@RMM
  overlandspinup = GetIntDefault("OverlandFlowSpinUp", 0);

  int itop, k1, io, io1, ovlnd_flag;           //DOK
  int ioo;         //@RMM

  double dtmp, dx, dy, dz, vol, vol2, ffx, ffy, ffz;          //@RMM
  double diff, coeff, x_coeff, y_coeff, z_coeff, updir, sep;         //@RMM
  double prod, prod_rt, prod_no, prod_up, prod_val, prod_lo;
  double prod_der, prod_rt_der, prod_no_der, prod_up_der;
  double west_temp, east_temp, north_temp, south_temp;
  double lower_temp, upper_temp, o_temp = 0.0;
  double sym_west_temp, sym_east_temp, sym_south_temp, sym_north_temp;
  double sym_lower_temp, sym_upper_temp;
  double lower_cond, upper_cond;

  //@RMM : terms for gravity/terrain
  double x_dir_g, y_dir_g, z_dir_g, del_x_slope, del_y_slope, x_dir_g_c, y_dir_g_c;

  BCStruct    *bc_struct;
  GrGeomSolid *gr_domain = ProblemDataGrDomain(problem_data);
  double      *bc_patch_values;
  double value, den_d, dend_d;
  int         *fdir;
  int ipatch, ival;

  CommHandle  *handle;
  VectorUpdateCommHandle  *vector_update_handle;

  // Determine if an overland flow boundary condition is being used.
  // If so will use the analytic Jacobian.
  if (public_xtra->type == not_set)
  {
    // Default to simple
    public_xtra->type = simple;

    BCPressureData   *bc_pressure_data
      = ProblemDataBCPressureData(problem_data);
    int num_patches = BCPressureDataNumPatches(bc_pressure_data);

    if (num_patches > 0)
    {
      int i;
      for (i = 0; i < num_patches; i++)
      {
        int type = BCPressureDataType(bc_pressure_data, i);
        switch (type)
        {
          case 7:
          {
            public_xtra->type = overland_flow;
          }
          break;
          case 10:
          {
            public_xtra->type = overland_flow;
          }
          break;
          case 11:
          {
            public_xtra->type = overland_flow;
          }
          break;
        }
      }
    }
  }

  /*-----------------------------------------------------------------------
   * Allocate temp vectors
   *-----------------------------------------------------------------------*/
  density_der = NewVectorType(grid, 1, 1, vector_cell_centered);
  saturation_der = NewVectorType(grid, 1, 1, vector_cell_centered);

  /*-----------------------------------------------------------------------
   * reuse the temp vectors for both saturation and rel_perm calculations.
   *-----------------------------------------------------------------------*/
  rel_perm = saturation;
  rel_perm_der = saturation_der;

  /* Pass pressure values to neighbors.  */
  vector_update_handle = InitVectorUpdate(pressure, VectorUpdateAll);
  FinalizeVectorUpdate(vector_update_handle);

/* Define grid for surface contribution */
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KW);
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KE);
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KN);
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KS);
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KWns);
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KEns);
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KNns);
  InlineNewVectorType(grid2d, 1, 1, vector_cell_centered, KSns);

  InlineInitVector(KW, 0.0);
  InlineInitVector(KE, 0.0);
  InlineInitVector(KN, 0.0);
  InlineInitVector(KS, 0.0);
  InlineInitVector(KWns, 0.0);
  InlineInitVector(KEns, 0.0);
  InlineInitVector(KNns, 0.0);
  InlineInitVector(KSns, 0.0);

  // SGS set this to 1 since the off/on behavior does not work in
  // parallel.
  ovlnd_flag = 1;  // determines whether or not to set up data structs for overland flow contribution


  /* Initialize matrix values to zero. */
  InlineInitMatrix(J, 0.0);
  InlineInitMatrix(JC, 0.0);

  /* Calculate time term contributions. */


  /*
    PFModuleInvokeType(PhaseDensityInvoke, density_module, (0, pressure, density, &dtmp, &dtmp,
    CALCFCN));
    PFModuleInvokeType(PhaseDensityInvoke, density_module, (0, pressure, density_der, &dtmp,
    &dtmp, CALCDER));
  */

  {
    GetModulePublicXtra(PhaseDensity, density_module, pdensity_xtra);

    int pd_phase = 0;
    Vector *phase_pressure = pressure;
    Vector *density_v = density;
    double *pressure_d = &dtmp;
    double *density_d = &dtmp;
    int fcn = CALCFCN;
    int sg;

#if 0
    // Fused loop version, since we know neither density is null and use the same phase type
    // as well as share the same coordinates
    /*
      Inline Fused: 0.000278978 seconds on average (1.47x speedup!)
      Inline Split: 0.000407163 seconds on average (1.028x speedup)
      Inline:       0.000409753 seconds on average (1.028x speedup)
      Function Ptr: 0.000411579 seconds on average

      This suggests that inlining is beneficial, but not nearly as much as fusing
      The overhead from function pointers seems minimal if not paired with other, similar
      looping structures nearby.
      EX: The Saturation vector function calls use GrGeomInLoop, meaning if THOSE were
      inlined, we could fuse them with the first GrGeomInLoop call of RJE.
      problem_saturation has significantly complex control flow however, so that
      may be difficult to do.
    */
    switch ((pdensity_xtra->type[pd_phase])) {
      case 0:
      {
        PhaseDensityType0 *dummy = (PhaseDensityType0*)(pdensity_xtra->data[pd_phase]);
        const double constant = (dummy->constant);
        Grid *pd_grid = VectorGrid(density_v);
        ForSubgridI(sg, GridSubgrids(pd_grid))
        {
          subgrid = GridSubgrid(pd_grid, sg);

          d_sub = VectorSubvector(density_v, sg);
          dd_sub = VectorSubvector(density_der, sg);

          ix = SubgridIX(subgrid) - 1;
          iy = SubgridIY(subgrid) - 1;
          iz = SubgridIZ(subgrid) - 1;

          nx = SubgridNX(subgrid) + 2;
          ny = SubgridNY(subgrid) + 2;
          nz = SubgridNZ(subgrid) + 2;

          nx_d = SubvectorNX(d_sub);
          ny_d = SubvectorNY(d_sub);
          nz_d = SubvectorNZ(d_sub);

          dp = SubvectorElt(d_sub, ix, iy, iz);
          ddp = SubvectorElt(dd_sub, ix, iy, iz);

          id = 0;
          BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                    id, nx_d, ny_d, nz_d, 1, 1, 1,
          {
            dp[id] = constant;
            ddp[id] = 0.0;
          });
        }
        break;
      }

      case 1:
      {
        double ref, comp;
        PhaseDensityType1 *dummy = (PhaseDensityType1*)(pdensity_xtra->data[pd_phase]);
        ref = (dummy->reference_density);
        comp = (dummy->compressibility_constant);

        Grid *pd_grid = VectorGrid(density_v);
        ForSubgridI(sg, GridSubgrids(pd_grid))
        {
          subgrid = GridSubgrid(pd_grid, sg);

          p_sub = VectorSubvector(phase_pressure, sg);
          d_sub = VectorSubvector(density_v, sg);
          dd_sub = VectorSubvector(density_der, sg);

          ix = SubgridIX(subgrid) - 1;
          iy = SubgridIY(subgrid) - 1;
          iz = SubgridIZ(subgrid) - 1;

          nx = SubgridNX(subgrid) + 2;
          ny = SubgridNY(subgrid) + 2;
          nz = SubgridNZ(subgrid) + 2;

          nx_p = SubvectorNX(p_sub);
          ny_p = SubvectorNY(p_sub);
          nz_p = SubvectorNZ(p_sub);

          nx_d = SubvectorNX(d_sub);
          ny_d = SubvectorNY(d_sub);
          nz_d = SubvectorNZ(d_sub);

          pp = SubvectorElt(p_sub, ix, iy, iz);
          dp = SubvectorElt(d_sub, ix, iy, iz);
          ddp = SubvectorElt(dd_sub, ix, iy, iz);

          ip = 0;
          id = 0;
          BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                    ip, nx_p, ny_p, nz_p, 1, 1, 1,
                    id, nx_d, ny_d, nz_d, 1, 1, 1,
          {
            dp[id] = ref * exp(pp[ip] * comp);
            ddp[id] = comp * ref * exp(pp[ip] * comp);
          });
        }
        break;
      }
    }
#endif
#if 0
    // Simpler version, since we know density_v is not null
    switch ((pdensity_xtra->type[pd_phase])) {
      case 0:
      {
        PhaseDensityType0 *dummy = (PhaseDensityType0*)(pdensity_xtra->data[pd_phase]);
        const double constant = (dummy->constant);
        Grid *pd_grid = VectorGrid(density_v);
        ForSubgridI(sg, GridSubgrids(pd_grid))
        {
          subgrid = GridSubgrid(pd_grid, sg);

          d_sub = VectorSubvector(density_v, sg);

          ix = SubgridIX(subgrid) - 1;
          iy = SubgridIY(subgrid) - 1;
          iz = SubgridIZ(subgrid) - 1;

          nx = SubgridNX(subgrid) + 2;
          ny = SubgridNY(subgrid) + 2;
          nz = SubgridNZ(subgrid) + 2;

          nx_d = SubvectorNX(d_sub);
          ny_d = SubvectorNY(d_sub);
          nz_d = SubvectorNZ(d_sub);

          dp = SubvectorElt(d_sub, ix, iy, iz);

          id = 0;
          BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                    id, nx_d, ny_d, nz_d, 1, 1, 1,
          {
            dp[id] = constant;
          });
        }

        density_v = density_der;
        pd_grid = VectorGrid(density_v);
        ForSubgridI(sg, GridSubgrids(pd_grid))
        {
          subgrid = GridSubgrid(pd_grid, sg);

          d_sub = VectorSubvector(density_v, sg);

          ix = SubgridIX(subgrid) - 1;
          iy = SubgridIY(subgrid) - 1;
          iz = SubgridIZ(subgrid) - 1;

          nx = SubgridNX(subgrid) + 2;
          ny = SubgridNY(subgrid) + 2;
          nz = SubgridNZ(subgrid) + 2;

          nx_d = SubvectorNX(d_sub);
          ny_d = SubvectorNY(d_sub);
          nz_d = SubvectorNZ(d_sub);

          dp = SubvectorElt(d_sub, ix, iy, iz);

          id = 0;
          BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                    id, nx_d, ny_d, nz_d, 1, 1, 1,
          {
            dp[id] = 0.0;
          });
        }
        break;
      }

      case 1:
      {
        double ref, comp;
        PhaseDensityType1 *dummy = (PhaseDensityType1*)(pdensity_xtra->data[pd_phase]);
        ref = (dummy->reference_density);
        comp = (dummy->compressibility_constant);

        Grid *pd_grid = VectorGrid(density_v);
        ForSubgridI(sg, GridSubgrids(pd_grid))
        {
          subgrid = GridSubgrid(pd_grid, sg);

          p_sub = VectorSubvector(phase_pressure, sg);
          d_sub = VectorSubvector(density_v, sg);

          ix = SubgridIX(subgrid) - 1;
          iy = SubgridIY(subgrid) - 1;
          iz = SubgridIZ(subgrid) - 1;

          nx = SubgridNX(subgrid) + 2;
          ny = SubgridNY(subgrid) + 2;
          nz = SubgridNZ(subgrid) + 2;

          nx_p = SubvectorNX(p_sub);
          ny_p = SubvectorNY(p_sub);
          nz_p = SubvectorNZ(p_sub);

          nx_d = SubvectorNX(d_sub);
          ny_d = SubvectorNY(d_sub);
          nz_d = SubvectorNZ(d_sub);

          pp = SubvectorElt(p_sub, ix, iy, iz);
          dp = SubvectorElt(d_sub, ix, iy, iz);

          ip = 0;
          id = 0;
          BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                    ip, nx_p, ny_p, nz_p, 1, 1, 1,
                    id, nx_d, ny_d, nz_d, 1, 1, 1,
          {
            dp[id] = ref * exp(pp[ip] * comp);
          });
        }

        density_v = density_der;
        pd_grid = VectorGrid(density_v);
        ForSubgridI(sg, GridSubgrids(pd_grid))
        {
          subgrid = GridSubgrid(pd_grid, sg);

          p_sub = VectorSubvector(phase_pressure, sg);
          d_sub = VectorSubvector(density_v, sg);

          ix = SubgridIX(subgrid) - 1;
          iy = SubgridIY(subgrid) - 1;
          iz = SubgridIZ(subgrid) - 1;

          nx = SubgridNX(subgrid) + 2;
          ny = SubgridNY(subgrid) + 2;
          nz = SubgridNZ(subgrid) + 2;

          nx_p = SubvectorNX(p_sub);
          ny_p = SubvectorNY(p_sub);
          nz_p = SubvectorNZ(p_sub);

          nx_d = SubvectorNX(d_sub);
          ny_d = SubvectorNY(d_sub);
          nz_d = SubvectorNZ(d_sub);

          pp = SubvectorElt(p_sub, ix, iy, iz);
          dp = SubvectorElt(d_sub, ix, iy, iz);

          ip = 0;
          id = 0;
          BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                    ip, nx_p, ny_p, nz_p, 1, 1, 1,
                    id, nx_d, ny_d, nz_d, 1, 1, 1,
          {
            dp[id] = comp * ref * exp(pp[ip] * comp);
          });
        }
        break;
      }
    }
#endif
#if 1
    switch ((pdensity_xtra->type[pd_phase])) {
      case 0:
      {
        double constant;
        PhaseDensityType0 *dummy = (PhaseDensityType0*)(pdensity_xtra->data[pd_phase]);
        constant = (dummy->constant);

        if (density_v != NULL)
        {
          Grid *pd_grid = VectorGrid(density_v);
          ForSubgridI(sg, GridSubgrids(pd_grid))
          {
            subgrid = GridSubgrid(pd_grid, sg);

            d_sub = VectorSubvector(density_v, sg);

            ix = SubgridIX(subgrid) - 1;
            iy = SubgridIY(subgrid) - 1;
            iz = SubgridIZ(subgrid) - 1;

            nx = SubgridNX(subgrid) + 2;
            ny = SubgridNY(subgrid) + 2;
            nz = SubgridNZ(subgrid) + 2;

            nx_d = SubvectorNX(d_sub);
            ny_d = SubvectorNY(d_sub);
            nz_d = SubvectorNZ(d_sub);

            dp = SubvectorElt(d_sub, ix, iy, iz);

            id = 0;
            if (fcn == CALCFCN)
            {
              BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = constant;
              });
            }
            else   /* fcn = CALCDER */
            {
              BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = 0.0;
              });
            }   /* End if fcn */
          }    /* End subgrid loop */
        }      /* End if density_v is not NULL */
        else
        {
          if (fcn == CALCFCN)
          {
            (*density_d) = constant;
          }
          else  /* fcn = CALCDER */
          {
            (*density_d) = 0.0;
          }
        }      /* End else */
        break;
      }        /* End case 0 */

      case 1:
      {
        double ref, comp;
        PhaseDensityType1 *dummy = (PhaseDensityType1*)(pdensity_xtra->data[pd_phase]);
        ref = (dummy->reference_density);
        comp = (dummy->compressibility_constant);

        if (density_v != NULL)
        {
          Grid *pd_grid = VectorGrid(density_v);
          ForSubgridI(sg, GridSubgrids(pd_grid))
          {
            subgrid = GridSubgrid(pd_grid, sg);

            p_sub = VectorSubvector(phase_pressure, sg);
            d_sub = VectorSubvector(density_v, sg);

            ix = SubgridIX(subgrid) - 1;
            iy = SubgridIY(subgrid) - 1;
            iz = SubgridIZ(subgrid) - 1;

            nx = SubgridNX(subgrid) + 2;
            ny = SubgridNY(subgrid) + 2;
            nz = SubgridNZ(subgrid) + 2;

            nx_p = SubvectorNX(p_sub);
            ny_p = SubvectorNY(p_sub);
            nz_p = SubvectorNZ(p_sub);

            nx_d = SubvectorNX(d_sub);
            ny_d = SubvectorNY(d_sub);
            nz_d = SubvectorNZ(d_sub);

            pp = SubvectorElt(p_sub, ix, iy, iz);
            dp = SubvectorElt(d_sub, ix, iy, iz);

            ip = 0;
            id = 0;

            if (fcn == CALCFCN)
            {
              BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                        ip, nx_p, ny_p, nz_p, 1, 1, 1,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = ref * exp(pp[ip] * comp);
              });
            }
            else          /* fcn = CALCDER */
            {
              BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                        ip, nx_p, ny_p, nz_p, 1, 1, 1,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = comp * ref * exp(pp[ip] * comp);
              });
            }
          }
        }
        else
        {
          if (fcn == CALCFCN)
          {
            (*density_d) = ref * exp((*pressure_d) * comp);
          }
          else
          {
            (*density_d) = comp * ref * exp((*pressure_d) * comp);
          }
        }

        break;
      }        /* End case 1 */
    }


    // Do it again for density_der
    density_v = density_der;
    fcn = CALCDER;
    switch ((pdensity_xtra->type[pd_phase])) {
      case 0:
      {
        double constant;
        PhaseDensityType0 *dummy = (PhaseDensityType0*)(pdensity_xtra->data[pd_phase]);
        constant = (dummy->constant);

        if (density_v != NULL)
        {
          Grid *pd_grid = VectorGrid(density_v);
          ForSubgridI(sg, GridSubgrids(pd_grid))
          {
            subgrid = GridSubgrid(pd_grid, sg);

            d_sub = VectorSubvector(density_v, sg);

            ix = SubgridIX(subgrid) - 1;
            iy = SubgridIY(subgrid) - 1;
            iz = SubgridIZ(subgrid) - 1;

            nx = SubgridNX(subgrid) + 2;
            ny = SubgridNY(subgrid) + 2;
            nz = SubgridNZ(subgrid) + 2;

            nx_d = SubvectorNX(d_sub);
            ny_d = SubvectorNY(d_sub);
            nz_d = SubvectorNZ(d_sub);

            dp = SubvectorElt(d_sub, ix, iy, iz);

            id = 0;
            if (fcn == CALCFCN)
            {
              BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = constant;
              });
            }
            else   /* fcn = CALCDER */
            {
              BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = 0.0;
              });
            }   /* End if fcn */
          }    /* End subgrid loop */
        }      /* End if density_v is not NULL */
        else
        {
          if (fcn == CALCFCN)
          {
            (*density_d) = constant;
          }
          else  /* fcn = CALCDER */
          {
            (*density_d) = 0.0;
          }
        }      /* End else */
        break;
      }        /* End case 0 */

      case 1:
      {
        double ref, comp;
        PhaseDensityType1 *dummy = (PhaseDensityType1*)(pdensity_xtra->data[pd_phase]);
        ref = (dummy->reference_density);
        comp = (dummy->compressibility_constant);

        if (density_v != NULL)
        {
          Grid *pd_grid = VectorGrid(density_v);
          ForSubgridI(sg, GridSubgrids(pd_grid))
          {
            subgrid = GridSubgrid(pd_grid, sg);

            p_sub = VectorSubvector(phase_pressure, sg);
            d_sub = VectorSubvector(density_v, sg);

            ix = SubgridIX(subgrid) - 1;
            iy = SubgridIY(subgrid) - 1;
            iz = SubgridIZ(subgrid) - 1;

            nx = SubgridNX(subgrid) + 2;
            ny = SubgridNY(subgrid) + 2;
            nz = SubgridNZ(subgrid) + 2;

            nx_p = SubvectorNX(p_sub);
            ny_p = SubvectorNY(p_sub);
            nz_p = SubvectorNZ(p_sub);

            nx_d = SubvectorNX(d_sub);
            ny_d = SubvectorNY(d_sub);
            nz_d = SubvectorNZ(d_sub);

            pp = SubvectorElt(p_sub, ix, iy, iz);
            dp = SubvectorElt(d_sub, ix, iy, iz);

            ip = 0;
            id = 0;

            if (fcn == CALCFCN)
            {
              BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                        ip, nx_p, ny_p, nz_p, 1, 1, 1,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = ref * exp(pp[ip] * comp);
              });
            }
            else          /* fcn = CALCDER */
            {
              BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                        ip, nx_p, ny_p, nz_p, 1, 1, 1,
                        id, nx_d, ny_d, nz_d, 1, 1, 1,
              {
                dp[id] = comp * ref * exp(pp[ip] * comp);
              });
            }
          }
        }
        else
        {
          if (fcn == CALCFCN)
          {
            (*density_d) = ref * exp((*pressure_d) * comp);
          }
          else
          {
            (*density_d) = comp * ref * exp((*pressure_d) * comp);
          }
        }

        break;
      }        /* End case 1 */
    }
#endif
  }

  // Deal with scope variable name conflicts
  {
    GetModulePublicXtra(Saturation, saturation_module, sat_xtra);
    Vector *phase_saturation;
    Vector *phase_pressure;
    Vector *phase_density;
    //double gravity; // Already in scope
    //ProblemData *problem_data; // Already in scope
    int fcn, sg;

    phase_saturation = saturation;
    phase_pressure = pressure;
    phase_density = density;
    fcn = CALCFCN;

    GrGeomSolid   *gr_solid, *gr_domain;

    Subvector     *ps_sub;
    Subvector     *pp_sub;
    Subvector     *pd_sub;
    Subvector     *satRF_sub;
    Subvector     *n_values_sub;
    Subvector     *alpha_values_sub;
    Subvector     *s_res_values_sub;
    Subvector     *s_sat_values_sub;

    double        *psdat, *ppdat, *pddat, *satRFdat;
    double        *n_values_dat, *alpha_values_dat;
    double        *s_res_values_dat, *s_sat_values_dat;

    int ix, iy, iz, r;
    int nx, ny, nz;

    int i, j, k, ips, ipp, ipd, ipRF;

    int n_index, alpha_index, s_res_index, s_sat_index;

    int            *region_indices, num_regions, ir;

    Grid *grid = VectorGrid(phase_saturation);
    SubgridArray *subgrids = GridSubgrids(grid);
#undef max
    InitVectorAll(phase_saturation, -FLT_MAX);

    SATURATION_MODULE;

    // Repeat for saturation_der
    fcn = CALCDER;
    phase_saturation = saturation_der;
    grid = VectorGrid(phase_saturation);
    subgrids = GridSubgrids(grid);
    InitVectorAll(phase_saturation, -FLT_MAX);

    SATURATION_MODULE;

  }


  /*
  PFModuleInvokeType(SaturationInvoke, saturation_module, (saturation, pressure,
                                                           density, gravity, problem_data,
                                                           CALCFCN));
  PFModuleInvokeType(SaturationInvoke, saturation_module, (saturation_der, pressure,
                                                           density, gravity, problem_data,
                                                           CALCDER));
  */

  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    J_sub = MatrixSubmatrix(J, is);
    cp = SubmatrixStencilData(J_sub, 0);

    p_sub = VectorSubvector(pressure, is);
    d_sub = VectorSubvector(density, is);
    s_sub = VectorSubvector(saturation, is);
    dd_sub = VectorSubvector(density_der, is);
    sd_sub = VectorSubvector(saturation_der, is);
    po_sub = VectorSubvector(porosity, is);
    ss_sub = VectorSubvector(sstorage, is);

    /* @RMM added to provide access to zmult */
    z_mult_sub = VectorSubvector(z_mult, is);
    /* @RMM added to provide variable dz */
    z_mult_dat = SubvectorData(z_mult_sub);
    /* @RMM added to provide access to x/y slopes */
    x_ssl_sub = VectorSubvector(x_ssl, is);
    y_ssl_sub = VectorSubvector(y_ssl, is);
    /* @RMM  added to provide slopes to terrain fns */
    x_ssl_dat = SubvectorData(x_ssl_sub);
    y_ssl_dat = SubvectorData(y_ssl_sub);

    /* RDF: assumes resolutions are the same in all 3 directions */
    r = SubgridRX(subgrid);

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);
    iz = SubgridIZ(subgrid);

    nx = SubgridNX(subgrid);
    ny = SubgridNY(subgrid);
    nz = SubgridNZ(subgrid);

    dx = SubgridDX(subgrid);
    dy = SubgridDY(subgrid);
    dz = SubgridDZ(subgrid);

    vol = dx * dy * dz;

    nx_v = SubvectorNX(d_sub);
    ny_v = SubvectorNY(d_sub);
    nz_v = SubvectorNZ(d_sub);

    nx_po = SubvectorNX(po_sub);
    ny_po = SubvectorNY(po_sub);
    nz_po = SubvectorNZ(po_sub);

    nx_m = SubmatrixNX(J_sub);
    ny_m = SubmatrixNY(J_sub);
    nz_m = SubmatrixNZ(J_sub);

    pp = SubvectorData(p_sub);     //pressure
    dp = SubvectorData(d_sub);     // density
    sp = SubvectorData(s_sub);     //saturation
    ddp = SubvectorData(dd_sub);    // density derivative: del-rho / del-press
    sdp = SubvectorData(sd_sub);    // saturation derivative: del-S / del-press
    pop = SubvectorData(po_sub);     // porosity
    ss = SubvectorData(ss_sub);     // sepcific storage

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      im = SubmatrixEltIndex(J_sub, i, j, k);
      ipo = SubvectorEltIndex(po_sub, i, j, k);
      iv = SubvectorEltIndex(d_sub, i, j, k);
      vol2 = vol * z_mult_dat[ipo];
      cp[im] += (sdp[iv] * dp[iv] + sp[iv] * ddp[iv])
                * pop[ipo] * vol2 + ss[iv] * vol2 * (sdp[iv] * dp[iv] * pp[iv] + sp[iv] * ddp[iv] * pp[iv] + sp[iv] * dp[iv]); //sk start
    });
  }    /* End subgrid loop */


  {
    GetModulePublicXtra(BCPressure, bc_pressure, public_xtra);
    GetModuleInstanceXtra(BCPressure, bc_pressure, instance_xtra);

    PFModule       *phase_density = (instance_xtra->phase_density);

    BCPressureData *bc_pressure_data = ProblemDataBCPressureData(problem_data);

    TimeCycleData  *time_cycle_data;

    int num_phases = (public_xtra->num_phases);

    Problem        *problem = (instance_xtra->problem);

    SubgridArray   *subgrids = GridSubgrids(grid);

    Subgrid        *subgrid;

    Vector      *z_mult = ProblemDataZmult(problem_data);
    Vector      *rsz = ProblemDataRealSpaceZ(problem_data);
    Subvector   *z_mult_sub;
    Subvector   *rsz_sub;
    double      *z_mult_dat;
    double      *rsz_dat;

    double       ***values;

    double         *patch_values;
    int patch_values_size;

    int            *fdir;

    int num_patches;
    int ipatch, is, i, j, k, ival, ips, phase;
    int cycle_number, interval_number;


    //       if (time == 10000.0) {
    //      printf("time: %f \n", time);
    //bc_struct = NULL; // Already in scope

    num_patches = BCPressureDataNumPatches(bc_pressure_data);

    if (num_patches > 0)
    {
      time_cycle_data = BCPressureDataTimeCycleData(bc_pressure_data);

      /*---------------------------------------------------------------------
       * Set up bc_struct with NULL values component
       *---------------------------------------------------------------------*/

      bc_struct = NewBCStruct(subgrids, gr_domain,
                              num_patches,
                              BCPressureDataPatchIndexes(bc_pressure_data),
                              BCPressureDataBCTypes(bc_pressure_data),
                              NULL);

      /*---------------------------------------------------------------------
       * Set up values component of bc_struct
       *---------------------------------------------------------------------*/

      values = ctalloc(double **, num_patches);
      BCStructValues(bc_struct) = values;

      for (ipatch = 0; ipatch < num_patches; ipatch++)
      {
        values[ipatch] = ctalloc(double *, SubgridArraySize(subgrids));

        cycle_number = BCPressureDataCycleNumber(bc_pressure_data, ipatch);
        interval_number = TimeCycleDataComputeIntervalNumber(
                                                             problem, time, time_cycle_data, cycle_number);

        switch (BCPressureDataType(bc_pressure_data, ipatch))
        {
          case DirEquilRefPatch:
          {
            /* Constant pressure value specified on a reference patch.
             * Calculate hydrostatic conditions along boundary patch for
             * elevations different from reference patch elevations.
             * Hydrostatic condition is:
             * grad p - rho g grad z = 0 */

            GeomSolid       *ref_solid;
            double z, dz2, dtmp;
            double offset, interface_press, interface_den;
            double ref_den, ref_press, nonlin_resid;
            double density_der, density, fcn_val;
            double height;
            double gravity = -ProblemGravity(problem);

            int ref_patch;
            int max_its = 10;
            int iterations;
            int ix, iy, iz, nx, ny, nz, r, iel;

            double         **elevations;

            GetBCPressureTypeStruct(DirEquilRefPatch, interval_data, bc_pressure_data,
                                    ipatch, interval_number);

            if (instance_xtra->elevations == NULL)
            {
              instance_xtra->elevations = ctalloc(double **, num_patches);
              instance_xtra->problem_data = problem_data;
              instance_xtra->grid = grid;
            }

            if (instance_xtra->elevations[ipatch] == NULL)
            {
              ref_solid = ProblemDataSolid(problem_data,
                                           DirEquilRefPatchRefSolid(interval_data));
              ref_patch = DirEquilRefPatchRefPatch(interval_data);

              /* Calculate elevations at (x,y) points on reference patch. */
              instance_xtra->elevations[ipatch] = CalcElevations(ref_solid, ref_patch, subgrids, problem_data);
            }

            elevations = instance_xtra->elevations[ipatch];

            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              z_mult_sub = VectorSubvector(z_mult, is);
              rsz_sub = VectorSubvector(rsz, is);
              z_mult_dat = SubvectorData(z_mult_sub);
              rsz_dat = SubvectorData(rsz_sub);

              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              ix = SubgridIX(subgrid);
              iy = SubgridIY(subgrid);
              iz = SubgridIZ(subgrid);

              nx = SubgridNX(subgrid);
              ny = SubgridNY(subgrid);
              nz = SubgridNZ(subgrid);

              /* RDF: assume resolution is the same in all 3 directions */
              r = SubgridRX(subgrid);

              dz2 = SubgridDZ(subgrid) * 0.5;

              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct,
                                ipatch, is,
              {
                ref_press = DirEquilRefPatchValue(interval_data);
                PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                   (0, NULL, NULL, &ref_press, &ref_den,
                                    CALCFCN));
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];
                iel = (i - ix) + (j - iy) * nx;
                fcn_val = 0.0;
                nonlin_resid = 1.0;
                iterations = -1;

                /* Solve a nonlinear problem for hydrostatic pressure
                 * at points on boundary patch given pressure on reference
                 * patch.  Note that the problem is only nonlinear if
                 * density depends on pressure.
                 *
                 * The nonlinear problem to solve is:
                 *   F(p) = 0
                 *   F(p) = P - P_ref
                 *          - 0.5*(rho(P) + rho(P_ref))*gravity*(z - z_ref)
                 *
                 * Newton's method is used to find a solution. */

                while ((nonlin_resid > 1.0E-6) && (iterations < max_its))
                {
                  if (iterations > -1)
                  {
                    PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                       (0, NULL, NULL, &patch_values[ival],
                                        &density_der, CALCDER));
                    dtmp = 1.0 - 0.5 * density_der * gravity
                           * (z - elevations[is][iel]);
                    patch_values[ival] = patch_values[ival] - fcn_val / dtmp;
                  }
                  else
                  {
                    patch_values[ival] = ref_press;
                  }
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                     (0, NULL, NULL, &patch_values[ival],
                                      &density, CALCFCN));

                  fcn_val = patch_values[ival] - ref_press
                            - 0.5 * (density + ref_den) * gravity
                            * (z - elevations[is][iel]);
                  nonlin_resid = fabs(fcn_val);

                  iterations++;
                }            /* End of while loop */


                /* Iterate over the phases and reset pressures according to
                 * hydrostatic conditions with appropriate densities.
                 * At each interface, we have hydrostatic conditions, so
                 *
                 * z_inter = (P_inter - P_ref) /
                 *            (0.5*(rho(P_inter)+rho(P_ref))*gravity
                 + z_ref
                 +
                 + Thus, the interface height and pressure are known
                 + and hydrostatic conditions can be determined for
                 + new phase.
                 +
                 + NOTE:  This only works for Pc = 0. */

                for (phase = 1; phase < num_phases; phase++)
                {
                  interface_press = DirEquilRefPatchValueAtInterface(
                                                                     interval_data, phase);
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                     (phase - 1, NULL, NULL, &interface_press,
                                      &interface_den, CALCFCN));
                  offset = (interface_press - ref_press)
                           / (0.5 * (interface_den + ref_den) * gravity);
                  ref_press = interface_press;
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                     (phase, NULL, NULL, &ref_press, &ref_den,
                                      CALCFCN));

                  /* Only reset pressure value if in another phase.
                   * The following "if" test determines whether this point
                   * is in another phase by checking if the computed
                   * pressure is less than the interface value.  This
                   * test ONLY works if the phases are distributed such
                   * that the lighter phases are above the heavier ones. */

                  if (patch_values[ival] < interface_press)
                  {
                    height = elevations[is][iel];
                    nonlin_resid = 1.0;
                    iterations = -1;
                    while ((nonlin_resid > 1.0E-6) && (iterations < max_its))
                    {
                      if (iterations > -1)
                      {
                        PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                           (phase, NULL, NULL, &patch_values[ival],
                                            &density_der, CALCDER));

                        dtmp = 1.0 - 0.5 * density_der * gravity
                               * (z - height);
                        patch_values[ival] = patch_values[ival]
                                             - fcn_val / dtmp;
                      }
                      else
                      {
                        height = height + offset;
                        patch_values[ival] = ref_press;
                      }

                      PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                         (phase, NULL, NULL,
                                          &patch_values[ival], &density,
                                          CALCFCN));

                      fcn_val = patch_values[ival] - ref_press
                                - 0.5 * (density + ref_den)
                                * gravity * (z - height);
                      nonlin_resid = fabs(fcn_val);

                      iterations++;
                    }              /* End of while loop */
                  }                /* End if above interface */
                }                  /* End phase loop */
              });                  /* End BCStructPatchLoop body */
            }                      /* End subgrid loop */


            break;
          } /* End DirEquilRefPatch */

          case DirEquilPLinear:
          {
            /* Piecewise linear pressure value specified on reference
             * patch.
             * Calculate hydrostatic conditions along patch for
             * elevations different from reference patch elevations.
             * Hydrostatic condition is:
             *             grad p - rho g grad z = 0 */

            int num_points;
            int ip;

            double x, y, z, dx2, dy2, dz2;
            double unitx, unity, line_min, line_length, xy, slope;

            double dtmp, offset, interface_press, interface_den;
            double ref_den, ref_press, nonlin_resid;
            double density_der, density, fcn_val;
            double height;
            double gravity = -ProblemGravity(problem);

            int max_its = 10;
            int iterations;

            GetBCPressureTypeStruct(DirEquilPLinear, interval_data, bc_pressure_data,
                                    ipatch, interval_number);

            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              z_mult_sub = VectorSubvector(z_mult, is);
              rsz_sub = VectorSubvector(rsz, is);
              z_mult_dat = SubvectorData(z_mult_sub);
              rsz_dat = SubvectorData(rsz_sub);

              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              dx2 = SubgridDX(subgrid) / 2.0;
              dy2 = SubgridDY(subgrid) / 2.0;
              dz2 = SubgridDZ(subgrid) / 2.0;

              /* compute unit direction vector for piecewise linear line */
              unitx = DirEquilPLinearXUpper(interval_data)
                      - DirEquilPLinearXLower(interval_data);
              unity = DirEquilPLinearYUpper(interval_data)
                      - DirEquilPLinearYLower(interval_data);
              line_length = sqrt(unitx * unitx + unity * unity);
              unitx /= line_length;
              unity /= line_length;
              line_min = DirEquilPLinearXLower(interval_data) * unitx
                         + DirEquilPLinearYLower(interval_data) * unity;

              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;
                y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2;
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];

                /* project center of BC face onto piecewise line */
                xy = (x * unitx + y * unity - line_min) / line_length;

                /* find two neighboring points */
                ip = 1;
                num_points = DirEquilPLinearNumPoints(interval_data);
                for (; ip < (num_points - 1); ip++)
                {
                  if (xy < DirEquilPLinearPoint(interval_data, ip))
                    break;
                }

                /* compute the slope */
                slope = ((DirEquilPLinearValue(interval_data, ip)
                          - DirEquilPLinearValue(interval_data, (ip - 1)))
                         / (DirEquilPLinearPoint(interval_data, ip)
                            - DirEquilPLinearPoint(interval_data, (ip - 1))));

                ref_press = DirEquilPLinearValue(interval_data, ip - 1)
                            + slope * (xy - DirEquilPLinearPoint(interval_data, ip - 1));
                PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                   (0, NULL, NULL, &ref_press, &ref_den,
                                    CALCFCN));
                fcn_val = 0.0;
                nonlin_resid = 1.0;
                iterations = -1;

                /* Solve a nonlinear problem for hydrostatic pressure
                 * at points on boundary patch given reference pressure.
                 * Note that the problem is only nonlinear if
                 * density depends on pressure.
                 *
                 * The nonlinear problem to solve is:
                 *   F(p) = 0
                 *   F(p) = P - P_ref
                 *          - 0.5*(rho(P) + rho(P_ref))*gravity*z
                 *
                 * Newton's method is used to find a solution. */

                while ((nonlin_resid > 1.0E-6) && (iterations < max_its))
                {
                  if (iterations > -1)
                  {
                    PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                       (0, NULL, NULL, &patch_values[ival],
                                        &density_der, CALCDER));
                    dtmp = 1.0 - 0.5 * density_der * gravity * z;
                    patch_values[ival] = patch_values[ival] - fcn_val / dtmp;
                  }
                  else
                  {
                    patch_values[ival] = ref_press;
                  }
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                     (0, NULL, NULL, &patch_values[ival],
                                      &density, CALCFCN));

                  fcn_val = patch_values[ival] - ref_press
                            - 0.5 * (density + ref_den) * gravity * z;
                  nonlin_resid = fabs(fcn_val);

                  iterations++;
                }            /* End of while loop */

                /* Iterate over the phases and reset pressures according to
                 * hydrostatic conditions with appropriate densities.
                 * At each interface, we have hydrostatic conditions, so
                 *
                 * z_inter = (P_inter - P_ref) /
                 *            (0.5*(rho(P_inter)+rho(P_ref))*gravity
                 + z_ref
                 +
                 + Thus, the interface height and pressure are known
                 + and hydrostatic conditions can be determined for
                 + new phase.
                 +
                 + NOTE:  This only works for Pc = 0. */

                for (phase = 1; phase < num_phases; phase++)
                {
                  interface_press = DirEquilPLinearValueAtInterface(
                                                                    interval_data, phase);
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                     (phase - 1, NULL, NULL, &interface_press,
                                      &interface_den, CALCFCN));
                  offset = (interface_press - ref_press)
                           / (0.5 * (interface_den + ref_den) * gravity);
                  ref_press = interface_press;
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                     (phase, NULL, NULL, &ref_press, &ref_den,
                                      CALCFCN));

                  /* Only reset pressure value if in another phase.
                   * The following "if" test determines whether this point
                   * is in another phase by checking if the computed
                   * pressure is less than the interface value.  This
                   * test ONLY works if the phases are distributed such
                   * that the lighter phases are above the heavier ones. */

                  if (patch_values[ival] < interface_press)
                  {
                    height = 0.0;
                    nonlin_resid = 1.0;
                    iterations = -1;
                    while ((nonlin_resid > 1.0E-6) && (iterations < max_its))
                    {
                      if (iterations > -1)
                      {
                        PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                           (phase, NULL, NULL, &patch_values[ival],
                                            &density_der, CALCDER));

                        dtmp = 1.0 - 0.5 * density_der * gravity * (z - height);
                        patch_values[ival] = patch_values[ival]
                                             - fcn_val / dtmp;
                      }
                      else
                      {
                        height = height + offset;
                        patch_values[ival] = ref_press;
                      }

                      PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                                         (phase, NULL, NULL,
                                          &patch_values[ival], &density,
                                          CALCFCN));

                      fcn_val = patch_values[ival] - ref_press
                                - 0.5 * (density + ref_den) * gravity
                                * (z - height);
                      nonlin_resid = fabs(fcn_val);

                      iterations++;
                    }              /* End of while loop */
                  }                /* End if above interface */
                }                  /* End phase loop */
              });                  /* End BCStructPatchLoop body */
            }
            break;
          } /* End DirEquilPLinear */

          case FluxConst:
          {
            /* Constant flux rate value on patch */
            double flux;

            GetBCPressureTypeStruct(FluxConst, interval_data, bc_pressure_data,
                                    ipatch, interval_number);

            flux = FluxConstValue(interval_data);
            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values[ival] = flux;
              });
            }       /* End subgrid loop */
            break;
          } /* End FluxConst */

          case FluxVolumetric:
          {
            /* Constant volumetric flux value on patch */
            double dx, dy, dz;
            double area, volumetric_flux;

            GetBCPressureTypeStruct(FluxVolumetric, interval_data, bc_pressure_data,
                                    ipatch, interval_number);


            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              z_mult_sub = VectorSubvector(z_mult, is);
              z_mult_dat = SubvectorData(z_mult_sub);

              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              dx = SubgridDX(subgrid);
              dy = SubgridDY(subgrid);
              dz = SubgridDZ(subgrid);

              area = 0.0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                /* primary direction x */
                if (fdir[0])
                {
                  area += dy * dz * z_mult_dat[ips];
                }
                /* primary direction y */
                else if (fdir[1])
                {
                  area += dx * dz * z_mult_dat[ips];
                }
                /* primary direction z */
                else if (fdir[2])
                {
                  area += dx * dy;
                }
              });

              if (area > 0.0)
              {
                volumetric_flux = FluxVolumetricValue(interval_data)
                                  / area;
                BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                {
                  patch_values[ival] = volumetric_flux;
                });
              }
            }           /* End subgrid loop */
            break;
          } /* End FluxVolumetric */

          case PressureFile:
          {
            /* Read input pressures from file (temporary).
             * This case assumes hydraulic head input conditions and
             * a constant density.  */
            Vector          *tmp_vector;
            Subvector       *subvector;
            char            *filename;
            double          *tmpp;
            int itmp;
            double z, dz2;
            double density, dtmp;

            double gravity = ProblemGravity(problem);

            /* Calculate density using dtmp as dummy argument. */
            dtmp = 0.0;
            PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                               (0, NULL, NULL, &dtmp, &density, CALCFCN));

            GetBCPressureTypeStruct(PressureFile, interval_data, bc_pressure_data,
                                    ipatch, interval_number);


            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              z_mult_sub = VectorSubvector(z_mult, is);
              rsz_sub = VectorSubvector(rsz, is);
              z_mult_dat = SubvectorData(z_mult_sub);
              rsz_dat = SubvectorData(rsz_sub);

              dz2 = SubgridDZ(subgrid) / 2.0;

              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              tmp_vector = NewVectorType(grid, 1, 0, vector_cell_centered);

              filename = PressureFileName(interval_data);
              ReadPFBinary(filename, tmp_vector);

              subvector = VectorSubvector(tmp_vector, is);

              tmpp = SubvectorData(subvector);
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];
                itmp = SubvectorEltIndex(subvector, i, j, k);

                patch_values[ival] = tmpp[itmp];     /*- density*gravity*z;*/
                /*last part taken out, very likely to be a bug)*/
              });

              FreeVector(tmp_vector);
            }             /* End subgrid loop */
            break;
          } /* End PressureFile */

          case FluxFile:
          {
            /* Read input fluxes from file (temporary) */
            Vector          *tmp_vector;
            Subvector       *subvector;
            char            *filename;
            double          *tmpp;
            int itmp;

            GetBCPressureTypeStruct(FluxFile, interval_data, bc_pressure_data,
                                    ipatch, interval_number);


            ForSubgridI(is, subgrids)
            {
              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              tmp_vector = NewVectorType(grid, 1, 0, vector_cell_centered);

              filename = FluxFileName(interval_data);
              ReadPFBinary(filename, tmp_vector);

              subvector = VectorSubvector(tmp_vector, is);

              tmpp = SubvectorData(subvector);
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                itmp = SubvectorEltIndex(subvector, i, j, k);

                patch_values[ival] = tmpp[itmp];
              });

              FreeVector(tmp_vector);
            }         /* End subgrid loop */
            break;
          } /* End FluxFile */

          case ExactSolution:
          {
            /* Calculate pressure based on pre-defined functions */
            double x, y, z, dx2, dy2, dz2;
            int fcn_type;

            GetBCPressureTypeStruct(ExactSolution, interval_data, bc_pressure_data,
                                    ipatch, interval_number);

            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              z_mult_sub = VectorSubvector(z_mult, is);
              rsz_sub = VectorSubvector(rsz, is);
              z_mult_dat = SubvectorData(z_mult_sub);
              rsz_dat = SubvectorData(rsz_sub);

              /* compute patch_values_size */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              dx2 = SubgridDX(subgrid) / 2.0;
              dy2 = SubgridDY(subgrid) / 2.0;
              dz2 = SubgridDZ(subgrid) / 2.0;

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              fcn_type = ExactSolutionFunctionType(interval_data);

              switch (fcn_type)
              {
                case 1:  /* p = x */
                {
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                  {
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;

                    patch_values[ival] = x;
                  });

                  break;
                }     /* End case 1 */

                case 2:  /* p = x + y + z */
                {
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                  {
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2;
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];
                    patch_values[ival] = x + y + z;
                  });

                  break;
                }     /* End case 2 */

                case 3:  /* p = x^3y^2 + sinxy + 1*/
                {
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                  {
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2;

                    patch_values[ival] = x * x * x * y * y + sin(x * y) + 1;
                  });
                  break;
                }     /* End case 3 */

                case 4:  /* p = x^3 y^4 + x^2 + sinxy cosy + 1 */
                {
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                  {
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2;
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];
                    patch_values[ival] = pow(x, 3) * pow(y, 4) + x * x + sin(x * y) * cos(y) + 1;
                  });
                  break;
                }     /* End case 4 */

                case 5:  /* p = xyzt +1 */
                {
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                  {
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2;
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];
                    patch_values[ival] = x * y * z * time + 1;
                  });
                  break;
                }     /* End case 5 */

                case 6:  /* p = xyzt +1 */
                {
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                  {
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2;
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];
                    patch_values[ival] = x * y * z * time + 1;
                  });
                  break;
                }     /* End case 5 */
              }       /* End switch */
            }         /* End subgrid loop */
            break;
          } /* End ExactSolution */

          case OverlandFlow:
          {
            /* Constant "rainfall" rate value on patch */
            double flux;

            GetBCPressureTypeStruct(FluxConst, interval_data, bc_pressure_data,
                                    ipatch, interval_number);


            flux = OverlandFlowValue(interval_data);
            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values[ival] = flux;
              });
            }       /* End subgrid loop */
            break;
          } /* End OverlandFlow */

          case OverlandFlowPFB:
          {
            /* Read input fluxes from file (overland) */
            Vector          *tmp_vector;
            Subvector       *subvector;
            //double          *data;
            char            *filename;
            double          *tmpp;
            int itmp;
            double dtmp;

            GetBCPressureTypeStruct(OverlandFlowPFB, interval_data, bc_pressure_data,
                                    ipatch, interval_number);

            ForSubgridI(is, subgrids)
            {
              /* compute patch_values_size (this isn't really needed yet) */
              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              tmp_vector = NewVectorType(grid, 1, 0, vector_cell_centered);
              //data = ctalloc(double, SizeOfVector(tmp_vector));
              //SetTempVectorData(tmp_vector, data);

              printf("reading overland file \n");
              filename = OverlandFlowPFBFileName(interval_data);
              ReadPFBinary(filename, tmp_vector);

              subvector = VectorSubvector(tmp_vector, is);

              tmpp = SubvectorData(subvector);
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                itmp = SubvectorEltIndex(subvector, i, j, k);

                patch_values[ival] = tmpp[itmp];
              });

              //tfree(VectorData(tmp_vector));
              FreeVector(tmp_vector);
            }              /* End subgrid loop */
            break;
          } /* End OverlandFlowPFB */

          case SeepageFace:
          {
            GetBCPressureTypeStruct(SeepageFace, interval_data, bc_pressure_data,
                                    ipatch, interval_number);
            double flux;

            flux = SeepageFaceValue(interval_data);
            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values[ival] = flux;
              });
            }

            break;
          } /* End SeepageFace */

          case OverlandKinematic:
          {
            GetBCPressureTypeStruct(OverlandKinematic, interval_data, bc_pressure_data,
                                    ipatch, interval_number);
            double flux;

            flux = OverlandKinematicValue(interval_data);
            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values[ival] = flux;
              });
            }

            break;
          } /* End OverlandKinematic */

          case OverlandDiffusive:
          {
            GetBCPressureTypeStruct(OverlandDiffusive, interval_data, bc_pressure_data,
                                    ipatch, interval_number);
            double flux;

            flux = OverlandDiffusiveValue(interval_data);
            ForSubgridI(is, subgrids)
            {
              subgrid = SubgridArraySubgrid(subgrids, is);

              patch_values_size = 0;
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values_size++;
              });

              patch_values = ctalloc(double, patch_values_size);
              values[ipatch][is] = patch_values;

              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                patch_values[ival] = flux;
              });
            }

            break;
          } /* End OverlandDiffusive */
        }
      }
    }
  }

  /*
  bc_struct = PFModuleInvokeType(BCPressureInvoke, bc_pressure,
                                 (problem_data, grid, gr_domain, time));
  */


  /* Get boundary pressure values for Dirichlet boundaries.   */
  /* These are needed for upstream weighting in mobilities - need boundary */
  /* values for rel perms and densities. */

  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    p_sub = VectorSubvector(pressure, is);

    nx_v = SubvectorNX(p_sub);
    ny_v = SubvectorNY(p_sub);
    nz_v = SubvectorNZ(p_sub);

    sy_v = nx_v;
    sz_v = ny_v * nx_v;

    pp = SubvectorData(p_sub);

    for (ipatch = 0; ipatch < BCStructNumPatches(bc_struct); ipatch++)
    {
      bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);

      switch (BCStructBCType(bc_struct, ipatch))
      {
        case DirichletBC:
        {
          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
          {
            ip = SubvectorEltIndex(p_sub, i, j, k);
            value = bc_patch_values[ival];
            pp[ip + fdir[0] * 1 + fdir[1] * sy_v + fdir[2] * sz_v] = value;
          });
          break;
        }
      }        /* End switch BCtype */
    }          /* End ipatch loop */
  }            /* End subgrid loop */


  // Calls problem_phase_rel_perm.c:PhaseRelPerm
  /* Calculate rel_perm and rel_perm_der */
  /* Arg list:
     Vector *phase_rel_perm
     Vector *phase_pressure
     Vector *phase_density
     double gravity
     ProblemData *problem_data
     int fcn
  */
  // These have static-inline VG table lookup calls in them, need to sort that out
  PFModuleInvokeType(PhaseRelPermInvoke, rel_perm_module,
                     (rel_perm, pressure, density, gravity, problem_data,
                      CALCFCN));

  PFModuleInvokeType(PhaseRelPermInvoke, rel_perm_module,
                     (rel_perm_der, pressure, density, gravity, problem_data,
                      CALCDER));

  /*
  {
    GetModulePublicXtra(PhaseRelPerm, rel_perm_module, public_xtra);
    Vector *phase_rel_perm = rel_perm;
    Vector *phase_pressure = pressure;
    Vector *phase_density = density;
    // double gravity = gravity; // Already in scope
    // ProblemData *problem_data = problem_data; // Already in scope
    int fcn = CALCFCN;

    PFModule      *this_module = ThisPFModule;
    PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

    Grid          *grid = VectorGrid(phase_rel_perm);

    GrGeomSolid   *gr_solid;

    Subvector     *pr_sub;
    Subvector     *pp_sub;
    Subvector     *pd_sub;
    Subvector     *n_values_sub;
    Subvector     *alpha_values_sub;

    double        *prdat, *ppdat, *pddat;
    double        *n_values_dat, *alpha_values_dat;

    SubgridArray  *subgrids = GridSubgrids(grid);

    Subgrid       *subgrid;

    int sg;

    int ix, iy, iz;
    int nx, ny, nz;

    int i, j, k, r, ipr, ipp, ipd;

    int n_index, alpha_index;

    int num_regions, *region_indices;
    int ir, *fdir;

    PHASE_REL_PERM_CALL;

    phase_rel_perm = rel_perm_der;
    fcn = CALCDER;

    PHASE_REL_PERM_CALL;
  }
  */

  /* Calculate contributions from second order derivatives and gravity */
  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);
    Subgrid* grid2d_subgrid = GridSubgrid(grid2d, is);
    int grid2d_iz = SubgridIZ(grid2d_subgrid);

    p_sub = VectorSubvector(pressure, is);
    d_sub = VectorSubvector(density, is);
    rp_sub = VectorSubvector(rel_perm, is);
    dd_sub = VectorSubvector(density_der, is);
    rpd_sub = VectorSubvector(rel_perm_der, is);
    permx_sub = VectorSubvector(permeability_x, is);
    permy_sub = VectorSubvector(permeability_y, is);
    permz_sub = VectorSubvector(permeability_z, is);
    J_sub = MatrixSubmatrix(J, is);

    /* @RMM added to provide access to x/y slopes */
    x_ssl_sub = VectorSubvector(x_ssl, is);
    y_ssl_sub = VectorSubvector(y_ssl, is);

    /* @RMM added to provide access to zmult */
    z_mult_sub = VectorSubvector(z_mult, is);
    /* @RMM added to provide variable dz */
    z_mult_dat = SubvectorData(z_mult_sub);

    r = SubgridRX(subgrid);

    ix = SubgridIX(subgrid) - 1;
    iy = SubgridIY(subgrid) - 1;
    iz = SubgridIZ(subgrid) - 1;

    nx = SubgridNX(subgrid) + 1;
    ny = SubgridNY(subgrid) + 1;
    nz = SubgridNZ(subgrid) + 1;

    dx = SubgridDX(subgrid);
    dy = SubgridDY(subgrid);
    dz = SubgridDZ(subgrid);

    ffx = dy * dz;
    ffy = dx * dz;
    ffz = dx * dy;

    nx_v = SubvectorNX(p_sub);
    ny_v = SubvectorNY(p_sub);
    nz_v = SubvectorNZ(p_sub);

    nx_m = SubmatrixNX(J_sub);
    ny_m = SubmatrixNY(J_sub);
    nz_m = SubmatrixNZ(J_sub);

    sy_v = nx_v;
    sz_v = ny_v * nx_v;
    sy_m = nx_m;
    sz_m = ny_m * nx_m;

    cp = SubmatrixStencilData(J_sub, 0);
    wp = SubmatrixStencilData(J_sub, 1);
    ep = SubmatrixStencilData(J_sub, 2);
    sop = SubmatrixStencilData(J_sub, 3);
    np = SubmatrixStencilData(J_sub, 4);
    lp = SubmatrixStencilData(J_sub, 5);
    up = SubmatrixStencilData(J_sub, 6);

    pp = SubvectorData(p_sub);
    dp = SubvectorData(d_sub);
    rpp = SubvectorData(rp_sub);
    ddp = SubvectorData(dd_sub);
    rpdp = SubvectorData(rpd_sub);
    permxp = SubvectorData(permx_sub);
    permyp = SubvectorData(permy_sub);
    permzp = SubvectorData(permz_sub);

    /* @RMM added to provide access FB values */
    FBx_sub = VectorSubvector(FBx, is);
    FBy_sub = VectorSubvector(FBy, is);
    FBz_sub = VectorSubvector(FBz, is);

    /* @RMM added to provide FB values */
    FBx_dat = SubvectorData(FBx_sub);
    FBy_dat = SubvectorData(FBy_sub);
    FBz_dat = SubvectorData(FBz_sub);

    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      ip = SubvectorEltIndex(p_sub, i, j, k);
      im = SubmatrixEltIndex(J_sub, i, j, k);
      ioo = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

      prod = rpp[ip] * dp[ip];
      prod_der = rpdp[ip] * dp[ip] + rpp[ip] * ddp[ip];

      prod_rt = rpp[ip + 1] * dp[ip + 1];
      prod_rt_der = rpdp[ip + 1] * dp[ip + 1] + rpp[ip + 1] * ddp[ip + 1];

      prod_no = rpp[ip + sy_v] * dp[ip + sy_v];
      prod_no_der = rpdp[ip + sy_v] * dp[ip + sy_v]
                    + rpp[ip + sy_v] * ddp[ip + sy_v];

      prod_up = rpp[ip + sz_v] * dp[ip + sz_v];
      prod_up_der = rpdp[ip + sz_v] * dp[ip + sz_v]
                    + rpp[ip + sz_v] * ddp[ip + sz_v];

      //@RMM  tfgupwind == 0 (default) should give original behavior
      // tfgupwind 1 should still use sine but upwind
      // tfgupwdin 2 just upwind
      switch (public_xtra->tfgupwind)
      {
        case 0:
          {
            // default formulation in Maxwell 2013
            x_dir_g = Mean(gravity * sin(atan(x_ssl_dat[ioo])), gravity * sin(atan(x_ssl_dat[ioo + 1])));
            x_dir_g_c = Mean(gravity * cos(atan(x_ssl_dat[ioo])), gravity * cos(atan(x_ssl_dat[ioo + 1])));
            y_dir_g = Mean(gravity * sin(atan(y_ssl_dat[ioo])), gravity * sin(atan(y_ssl_dat[ioo + sy_v])));
            y_dir_g_c = Mean(gravity * cos(atan(y_ssl_dat[ioo])), gravity * cos(atan(y_ssl_dat[ioo + sy_v])));
            break;
          }

        case 1:
          {
            // direct upwinding, no averaging with sines
            x_dir_g = gravity * sin(atan(x_ssl_dat[ioo]));
            x_dir_g_c = gravity * cos(atan(x_ssl_dat[ioo]));
            y_dir_g = gravity * sin(atan(y_ssl_dat[ioo]));
            y_dir_g_c = gravity * cos(atan(y_ssl_dat[ioo]));
            break;
          }

        case 2:
          {
            // direct upwinding, no averaging no sines
            x_dir_g = x_ssl_dat[ioo];
            x_dir_g_c = 1.0;
            y_dir_g = y_ssl_dat[ioo];
            y_dir_g_c = 1.0;
            break;
          }
      }


      /* diff >= 0 implies flow goes left to right */
      diff = pp[ip] - pp[ip + 1];
      updir = (diff / dx) * x_dir_g_c - x_dir_g;

      /* multiply X_coeff by FB in x */
      x_coeff = FBx_dat[ip] * dt * ffx * (1.0 / dx) * z_mult_dat[ip]
                * PMean(pp[ip], pp[ip + 1], permxp[ip], permxp[ip + 1])
                / viscosity;


      sym_west_temp = (-x_coeff
                       * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM TFG contributions, sym


      west_temp = (-x_coeff * diff
                   * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g_c
                  + sym_west_temp;

      west_temp += (x_coeff * dx * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g; //@RMM TFG contributions, non sym

      sym_east_temp = (-x_coeff
                       * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM added sym TFG contributions

      east_temp = (x_coeff * diff
                   * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g_c
                  + sym_east_temp;

      east_temp += -(x_coeff * dx * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g; //@RMM  TFG contributions non sym

      /* diff >= 0 implies flow goes south to north */
      diff = pp[ip] - pp[ip + sy_v];
      updir = (diff / dy) * y_dir_g_c - y_dir_g;


      /* multiply y_coeff by FB in y */
      y_coeff = FBy_dat[ip] * dt * ffy * (1.0 / dy) * z_mult_dat[ip]
                * PMean(pp[ip], pp[ip + sy_v], permyp[ip], permyp[ip + sy_v])
                / viscosity;

      sym_south_temp = -y_coeff
                       * RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM TFG contributions, SYMM

      south_temp = -y_coeff * diff
                   * RPMean(updir, 0.0, prod_der, 0.0) * y_dir_g_c
                   + sym_south_temp;

      south_temp += (y_coeff * dy * RPMean(updir, 0.0, prod_der, 0.0)) * y_dir_g; //@RMM TFG contributions, non sym


      sym_north_temp = y_coeff
                       * -RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM  TFG contributions non SYMM

      north_temp = y_coeff * diff
                   * RPMean(updir, 0.0, 0.0,
                            prod_no_der) * y_dir_g_c
                   + sym_north_temp;

      north_temp += -(y_coeff * dy * RPMean(updir, 0.0, 0.0, prod_no_der)) * y_dir_g; //@RMM  TFG contributions non sym

      sep = (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]));
      /* diff >= 0 implies flow goes lower to upper */


      lower_cond = pp[ip] / sep - (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip + sz_v])) * dp[ip] * gravity;

      upper_cond = pp[ip + sz_v] / sep + (z_mult_dat[ip + sz_v] / (z_mult_dat[ip] + z_mult_dat[ip + sz_v])) * dp[ip + sz_v] * gravity;


      diff = lower_cond - upper_cond;

      /* multiply z_coeff by FB in z */
      z_coeff = FBz_dat[ip] * dt * ffz
                * PMeanDZ(permzp[ip], permzp[ip + sz_v], z_mult_dat[ip], z_mult_dat[ip + sz_v])
                / viscosity;

      sym_lower_temp = -z_coeff * (1.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))
                       * RPMean(lower_cond, upper_cond, prod,
                                prod_up);

      lower_temp = -z_coeff
                   * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)
                      + (-gravity * 0.5 * dz * (Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])) * ddp[ip]
                         * RPMean(lower_cond, upper_cond, prod,
                                  prod_up)))
                   + sym_lower_temp;

      sym_upper_temp = z_coeff * (1.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))
                       * -RPMean(lower_cond, upper_cond, prod,
                                 prod_up);

      upper_temp = z_coeff
                   * (diff * RPMean(lower_cond, upper_cond, 0.0,
                                    prod_up_der)
                      + (-gravity * 0.5 * dz * (Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])) * ddp[ip + sz_v]
                         * RPMean(lower_cond, upper_cond, prod,
                                  prod_up)))
                   + sym_upper_temp;



      cp[im] -= west_temp + south_temp + lower_temp;
      cp[im + 1] -= east_temp;
      cp[im + sy_m] -= north_temp;
      cp[im + sz_m] -= upper_temp;

      if (!symm_part)
      {
        ep[im] += east_temp;
        np[im] += north_temp;
        up[im] += upper_temp;

        wp[im + 1] += west_temp;
        sop[im + sy_m] += south_temp;
        lp[im + sz_m] += lower_temp;
      }
      else     /* Symmetric matrix: just update upper coeffs */
      {
        ep[im] += sym_east_temp;
        np[im] += sym_north_temp;
        up[im] += sym_upper_temp;
      }
    });
  }  //

  /*  Calculate correction for boundary conditions */

  if (symm_part)
  {
    /*  For symmetric part only, we first adjust coefficients of normal */
    /*  direction boundary pressure by adding in the nonsymmetric part. */
    /*  The entire coefficicent will be subtracted from the diagonal    */
    /*  and set to zero in the subsequent section - no matter what type */
    /*  of BC is involved.  Without this correction, only the symmetric */
    /*  part would be removed, incorrectly leaving the nonsymmetric     */
    /*  contribution on the diagonal.                                   */

    ForSubgridI(is, GridSubgrids(grid))
    {
      subgrid = GridSubgrid(grid, is);

      p_sub = VectorSubvector(pressure, is);
      dd_sub = VectorSubvector(density_der, is);
      rpd_sub = VectorSubvector(rel_perm_der, is);
      d_sub = VectorSubvector(density, is);
      rp_sub = VectorSubvector(rel_perm, is);
      permx_sub = VectorSubvector(permeability_x, is);
      permy_sub = VectorSubvector(permeability_y, is);
      permz_sub = VectorSubvector(permeability_z, is);
      J_sub = MatrixSubmatrix(J, is);

      dx = SubgridDX(subgrid);
      dy = SubgridDY(subgrid);
      dz = SubgridDZ(subgrid);

      ffx = dy * dz;
      ffy = dx * dz;
      ffz = dx * dy;

      nx_v = SubvectorNX(p_sub);
      ny_v = SubvectorNY(p_sub);
      nz_v = SubvectorNZ(p_sub);

      sy_v = nx_v;
      sz_v = ny_v * nx_v;
      /* @RMM added to provide access to zmult */
      z_mult_sub = VectorSubvector(z_mult, is);
      /* @RMM added to provide variable dz */
      z_mult_dat = SubvectorData(z_mult_sub);

      cp = SubmatrixStencilData(J_sub, 0);
      wp = SubmatrixStencilData(J_sub, 1);
      ep = SubmatrixStencilData(J_sub, 2);
      sop = SubmatrixStencilData(J_sub, 3);
      np = SubmatrixStencilData(J_sub, 4);
      lp = SubmatrixStencilData(J_sub, 5);
      up = SubmatrixStencilData(J_sub, 6);

      pp = SubvectorData(p_sub);
      ddp = SubvectorData(dd_sub);
      rpdp = SubvectorData(rpd_sub);
      dp = SubvectorData(d_sub);
      rpp = SubvectorData(rp_sub);
      permxp = SubvectorData(permx_sub);
      permyp = SubvectorData(permy_sub);
      permzp = SubvectorData(permz_sub);

      for (ipatch = 0; ipatch < BCStructNumPatches(bc_struct); ipatch++)
      {
        BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
        {
          ip = SubvectorEltIndex(p_sub, i, j, k);
          im = SubmatrixEltIndex(J_sub, i, j, k);

          // SGS added this as prod was not being set to anything. Check with carol.
          prod = rpp[ip] * dp[ip];

          if (fdir[0])
          {
            switch (fdir[0])
            {
              case -1:
                {
                  diff = pp[ip - 1] - pp[ip];
                  prod_der = rpdp[ip - 1] * dp[ip - 1] + rpp[ip - 1] * ddp[ip - 1];
                  coeff = dt * z_mult_dat[ip] * ffx * (1.0 / dx)
                          * PMean(pp[ip - 1], pp[ip], permxp[ip - 1], permxp[ip])
                          / viscosity;
                  wp[im] = -coeff * diff
                           * RPMean(pp[ip - 1], pp[ip], prod_der, 0.0);
                  break;
                }

              case 1:
                {
                  diff = pp[ip] - pp[ip + 1];
                  prod_der = rpdp[ip + 1] * dp[ip + 1] + rpp[ip + 1] * ddp[ip + 1];
                  coeff = dt * z_mult_dat[ip] * ffx * (1.0 / dx)
                          * PMean(pp[ip], pp[ip + 1], permxp[ip], permxp[ip + 1])
                          / viscosity;
                  ep[im] = coeff * diff
                           * RPMean(pp[ip], pp[ip + 1], 0.0, prod_der);
                  break;
                }
            }         /* End switch on fdir[0] */
          }           /* End if (fdir[0]) */

          else if (fdir[1])
          {
            switch (fdir[1])
            {
              case -1:
                {
                  diff = pp[ip - sy_v] - pp[ip];
                  prod_der = rpdp[ip - sy_v] * dp[ip - sy_v]
                             + rpp[ip - sy_v] * ddp[ip - sy_v];
                  coeff = dt * z_mult_dat[ip] * ffy * (1.0 / dy)
                          * PMean(pp[ip - sy_v], pp[ip],
                                  permyp[ip - sy_v], permyp[ip])
                          / viscosity;
                  sop[im] = -coeff * diff
                            * RPMean(pp[ip - sy_v], pp[ip], prod_der, 0.0);

                  break;
                }

              case 1:
                {
                  diff = pp[ip] - pp[ip + sy_v];
                  prod_der = rpdp[ip + sy_v] * dp[ip + sy_v]
                             + rpp[ip + sy_v] * ddp[ip + sy_v];
                  coeff = dt * z_mult_dat[ip] * ffy * (1.0 / dy)
                          * PMean(pp[ip], pp[ip + sy_v],
                                  permyp[ip], permyp[ip + sy_v])
                          / viscosity;
                  np[im] = -coeff * diff
                           * RPMean(pp[ip], pp[ip + sy_v], 0.0, prod_der);
                  break;
                }
            }         /* End switch on fdir[1] */
          }           /* End if (fdir[1]) */

          else if (fdir[2])
          {
            switch (fdir[2])
            {
              case -1:
                {
                  lower_cond = (pp[ip - sz_v])
                               - 0.5 * dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v])
                               * dp[ip - sz_v] * gravity;
                  upper_cond = (pp[ip]) + 0.5 * dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v])
                               * dp[ip] * gravity;
                  diff = lower_cond - upper_cond;
                  prod_der = rpdp[ip - sz_v] * dp[ip - sz_v]
                             + rpp[ip - sz_v] * ddp[ip - sz_v];
                  prod_lo = rpp[ip - sz_v] * dp[ip - sz_v];
                  coeff = dt * ffz * (1.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v])))
                          * PMeanDZ(permzp[ip - sz_v], permzp[ip],
                                    z_mult_dat[ip - sz_v], z_mult_dat[ip])
                          / viscosity;
                  lp[im] = -coeff *
                           (diff * RPMean(lower_cond, upper_cond,
                                          prod_der, 0.0)
                            - gravity * 0.5 * dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v]) * ddp[ip]
                            * RPMean(lower_cond, upper_cond, prod_lo, prod));

                  break;
                }

              case 1:
                {
                  lower_cond = (pp[ip]) - 0.5 * dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * dp[ip] * gravity;
                  upper_cond = (pp[ip + sz_v])
                               + 0.5 * dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * dp[ip + sz_v] * gravity;
                  diff = lower_cond - upper_cond;
                  prod_der = rpdp[ip + sz_v] * dp[ip + sz_v]
                             + rpp[ip + sz_v] * ddp[ip + sz_v];
                  prod_up = rpp[ip + sz_v] * dp[ip + sz_v];
                  coeff = dt * ffz * (1.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))
                          * PMeanDZ(permzp[ip], permzp[ip + sz_v],
                                    z_mult_dat[ip], z_mult_dat[ip + sz_v])
                          / viscosity;
                  up[im] = -coeff *
                           (diff * RPMean(lower_cond, upper_cond,
                                          0.0, prod_der)
                            - gravity * 0.5 * dz * (Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])) * ddp[ip]
                            * RPMean(lower_cond, upper_cond, prod, prod_up));

                  break;
                }
            }         /* End switch on fdir[2] */
          }        /* End if (fdir[2]) */
        });       /* End Patch Loop */
      }           /* End ipatch loop */
    }             /* End subgrid loop */
  }                  /* End if symm_part */

  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    p_sub = VectorSubvector(pressure, is);
    s_sub = VectorSubvector(saturation, is);
    dd_sub = VectorSubvector(density_der, is);
    rpd_sub = VectorSubvector(rel_perm_der, is);
    d_sub = VectorSubvector(density, is);
    rp_sub = VectorSubvector(rel_perm, is);
    permx_sub = VectorSubvector(permeability_x, is);
    permy_sub = VectorSubvector(permeability_y, is);
    permz_sub = VectorSubvector(permeability_z, is);
    J_sub = MatrixSubmatrix(J, is);

    /* overland flow - DOK */
    kw_sub = VectorSubvector(KW, is);
    ke_sub = VectorSubvector(KE, is);
    kn_sub = VectorSubvector(KN, is);
    ks_sub = VectorSubvector(KS, is);
    kwns_sub = VectorSubvector(KWns, is);
    kens_sub = VectorSubvector(KEns, is);
    knns_sub = VectorSubvector(KNns, is);
    ksns_sub = VectorSubvector(KSns, is);

    dx = SubgridDX(subgrid);
    dy = SubgridDY(subgrid);
    dz = SubgridDZ(subgrid);

    /* @RMM added to provide access to zmult */
    z_mult_sub = VectorSubvector(z_mult, is);
    /* @RMM added to provide variable dz */
    z_mult_dat = SubvectorData(z_mult_sub);

    vol = dx * dy * dz;

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);

    ffx = dy * dz;
    ffy = dx * dz;
    ffz = dx * dy;

    nx_v = SubvectorNX(p_sub);
    ny_v = SubvectorNY(p_sub);
    nz_v = SubvectorNZ(p_sub);

    sy_v = nx_v;
    sz_v = ny_v * nx_v;

    cp = SubmatrixStencilData(J_sub, 0);
    wp = SubmatrixStencilData(J_sub, 1);
    ep = SubmatrixStencilData(J_sub, 2);
    sop = SubmatrixStencilData(J_sub, 3);
    np = SubmatrixStencilData(J_sub, 4);
    lp = SubmatrixStencilData(J_sub, 5);
    up = SubmatrixStencilData(J_sub, 6);

    /* overland flow contribution */
    kw_der = SubvectorData(kw_sub);
    ke_der = SubvectorData(ke_sub);
    kn_der = SubvectorData(kn_sub);
    ks_der = SubvectorData(ks_sub);
    kwns_der = SubvectorData(kwns_sub);
    kens_der = SubvectorData(kens_sub);
    knns_der = SubvectorData(knns_sub);
    ksns_der = SubvectorData(ksns_sub);

    pp = SubvectorData(p_sub);
    sp = SubvectorData(s_sub);
    ddp = SubvectorData(dd_sub);
    rpdp = SubvectorData(rpd_sub);
    dp = SubvectorData(d_sub);
    rpp = SubvectorData(rp_sub);
    permxp = SubvectorData(permx_sub);
    permyp = SubvectorData(permy_sub);
    permzp = SubvectorData(permz_sub);

    for (ipatch = 0; ipatch < BCStructNumPatches(bc_struct); ipatch++)
    {
      bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);

      switch (BCStructBCType(bc_struct, ipatch))
      {
        case DirichletBC:
        {
          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
          {
            ip = SubvectorEltIndex(p_sub, i, j, k);

            value = bc_patch_values[ival];

            // Arg names: phase, phase_pressure, density_v, pressure_d, density_d, fcn
            /*
            PFModuleInvokeType(PhaseDensityInvoke, density_module,
                               (0, NULL, NULL, &value, &den_d, CALCFCN));
            PFModuleInvokeType(PhaseDensityInvoke, density_module,
                               (0, NULL, NULL, &value, &dend_d, CALCDER));
            */
            // Variables redeclared as they appear in the above function calls

            {
              GetModulePublicXtra(PhaseDensity, density_module, public_xtra);

              Vector *phase_pressure = NULL;
              int phase = 0;
              Vector *density_v = NULL;
              double *pressure_d = &value;
              double *density_d = &den_d;
              int fcn = CALCFCN;

              double ref;
              double comp;

              Grid          *grid;
              Subvector     *p_sub;
              Subvector     *d_sub;
              double        *pp;
              double        *dp;
              Subgrid       *subgrid;
              int sg;
              int ix;
              int iy;
              int iz;
              int nx;
              int ny;
              int nz;
              int nx_p;
              int ny_p;
              int nz_p;
              int nx_d;
              int ny_d;
              int nz_d;
              int ip;
              int id;

              switch ((public_xtra->type[phase]))
              {
                case 0:
                {
                  double constant;
                  GetDummyType(PhaseDensity, 0, (public_xtra->data[phase]), dummy0);

                  constant = (dummy0->constant);

                  if (density_v != NULL)
                  {
                    grid = VectorGrid(density_v);
                    ForSubgridI(sg, GridSubgrids(grid))
                    {
                      subgrid = GridSubgrid(grid, sg);

                      d_sub = VectorSubvector(density_v, sg);

                      ix = SubgridIX(subgrid) - 1;
                      iy = SubgridIY(subgrid) - 1;
                      iz = SubgridIZ(subgrid) - 1;

                      nx = SubgridNX(subgrid) + 2;
                      ny = SubgridNY(subgrid) + 2;
                      nz = SubgridNZ(subgrid) + 2;

                      nx_d = SubvectorNX(d_sub);
                      ny_d = SubvectorNY(d_sub);
                      nz_d = SubvectorNZ(d_sub);

                      dp = SubvectorElt(d_sub, ix, iy, iz);

                      id = 0;
                      if (fcn == CALCFCN)
                      {
                        BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = constant;
                        });
                      }
                      else   /* fcn = CALCDER */
                      {
                        BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = 0.0;
                        });
                      }   /* End if fcn */
                    }    /* End subgrid loop */
                  }      /* End if density_v is not NULL */
                  else
                  {
                    if (fcn == CALCFCN)
                    {
                      (*density_d) = constant;
                    }
                    else  /* fcn = CALCDER */
                    {
                      (*density_d) = 0.0;
                    }
                  }      /* End else */
                  break;
                }        /* End case 0 */

                case 1:
                {
                  GetDummyType(PhaseDensity, 1, (public_xtra->data[phase]), dummy1);
                  ref = (dummy1->reference_density);
                  comp = (dummy1->compressibility_constant);

                  if (density_v != NULL)
                  {
                    grid = VectorGrid(density_v);
                    ForSubgridI(sg, GridSubgrids(grid))
                    {
                      subgrid = GridSubgrid(grid, sg);

                      p_sub = VectorSubvector(phase_pressure, sg);
                      d_sub = VectorSubvector(density_v, sg);

                      ix = SubgridIX(subgrid) - 1;
                      iy = SubgridIY(subgrid) - 1;
                      iz = SubgridIZ(subgrid) - 1;

                      nx = SubgridNX(subgrid) + 2;
                      ny = SubgridNY(subgrid) + 2;
                      nz = SubgridNZ(subgrid) + 2;

                      nx_p = SubvectorNX(p_sub);
                      ny_p = SubvectorNY(p_sub);
                      nz_p = SubvectorNZ(p_sub);

                      nx_d = SubvectorNX(d_sub);
                      ny_d = SubvectorNY(d_sub);
                      nz_d = SubvectorNZ(d_sub);

                      pp = SubvectorElt(p_sub, ix, iy, iz);
                      dp = SubvectorElt(d_sub, ix, iy, iz);

                      ip = 0;
                      id = 0;

                      if (fcn == CALCFCN)
                      {
                        BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                                  ip, nx_p, ny_p, nz_p, 1, 1, 1,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = ref * exp(pp[ip] * comp);
                        });
                      }
                      else          /* fcn = CALCDER */
                      {
                        BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                                  ip, nx_p, ny_p, nz_p, 1, 1, 1,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = comp * ref * exp(pp[ip] * comp);
                        });
                      }
                    }
                  }
                  else
                  {
                    if (fcn == CALCFCN)
                    {
                      (*density_d) = ref * exp((*pressure_d) * comp);
                    }
                    else
                    {
                      (*density_d) = comp * ref * exp((*pressure_d) * comp);
                    }
                  }

                  break;
                }        /* End case 1 */
              }

              fcn = CALCDER;
              switch ((public_xtra->type[phase]))
              {
                case 0:
                {
                  double constant;
                  GetDummyType(PhaseDensity, 0, (public_xtra->data[phase]), dummy0);
                  constant = (dummy0->constant);

                  if (density_v != NULL)
                  {
                    grid = VectorGrid(density_v);
                    ForSubgridI(sg, GridSubgrids(grid))
                    {
                      subgrid = GridSubgrid(grid, sg);

                      d_sub = VectorSubvector(density_v, sg);

                      ix = SubgridIX(subgrid) - 1;
                      iy = SubgridIY(subgrid) - 1;
                      iz = SubgridIZ(subgrid) - 1;

                      nx = SubgridNX(subgrid) + 2;
                      ny = SubgridNY(subgrid) + 2;
                      nz = SubgridNZ(subgrid) + 2;

                      nx_d = SubvectorNX(d_sub);
                      ny_d = SubvectorNY(d_sub);
                      nz_d = SubvectorNZ(d_sub);

                      dp = SubvectorElt(d_sub, ix, iy, iz);

                      id = 0;
                      if (fcn == CALCFCN)
                      {
                        BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = constant;
                        });
                      }
                      else   /* fcn = CALCDER */
                      {
                        BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = 0.0;
                        });
                      }   /* End if fcn */
                    }    /* End subgrid loop */
                  }      /* End if density_v is not NULL */
                  else
                  {
                    if (fcn == CALCFCN)
                    {
                      (*density_d) = constant;
                    }
                    else  /* fcn = CALCDER */
                    {
                      (*density_d) = 0.0;
                    }
                  }      /* End else */
                  break;
                }        /* End case 0 */

                case 1:
                {
                  GetDummyType(PhaseDensity, 1, (public_xtra->data[phase]), dummy1);
                  ref = (dummy1->reference_density);
                  comp = (dummy1->compressibility_constant);

                  if (density_v != NULL)
                  {
                    grid = VectorGrid(density_v);
                    ForSubgridI(sg, GridSubgrids(grid))
                    {
                      subgrid = GridSubgrid(grid, sg);

                      p_sub = VectorSubvector(phase_pressure, sg);
                      d_sub = VectorSubvector(density_v, sg);

                      ix = SubgridIX(subgrid) - 1;
                      iy = SubgridIY(subgrid) - 1;
                      iz = SubgridIZ(subgrid) - 1;

                      nx = SubgridNX(subgrid) + 2;
                      ny = SubgridNY(subgrid) + 2;
                      nz = SubgridNZ(subgrid) + 2;

                      nx_p = SubvectorNX(p_sub);
                      ny_p = SubvectorNY(p_sub);
                      nz_p = SubvectorNZ(p_sub);

                      nx_d = SubvectorNX(d_sub);
                      ny_d = SubvectorNY(d_sub);
                      nz_d = SubvectorNZ(d_sub);

                      pp = SubvectorElt(p_sub, ix, iy, iz);
                      dp = SubvectorElt(d_sub, ix, iy, iz);

                      ip = 0;
                      id = 0;

                      if (fcn == CALCFCN)
                      {
                        BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                                  ip, nx_p, ny_p, nz_p, 1, 1, 1,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = ref * exp(pp[ip] * comp);
                        });
                      }
                      else          /* fcn = CALCDER */
                      {
                        BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
                                  ip, nx_p, ny_p, nz_p, 1, 1, 1,
                                  id, nx_d, ny_d, nz_d, 1, 1, 1,
                        {
                          dp[id] = comp * ref * exp(pp[ip] * comp);
                        });
                      }
                    }
                  }
                  else
                  {
                    if (fcn == CALCFCN)
                    {
                      (*density_d) = ref * exp((*pressure_d) * comp);
                    }
                    else
                    {
                      (*density_d) = comp * ref * exp((*pressure_d) * comp);
                    }
                  }

                  break;
                }        /* End case 1 */
              }
            }


            ip = SubvectorEltIndex(p_sub, i, j, k);
            im = SubmatrixEltIndex(J_sub, i, j, k);

            prod = rpp[ip] * dp[ip];
            prod_der = rpdp[ip] * dp[ip] + rpp[ip] * ddp[ip];

            if (fdir[0])
            {
              coeff = dt * ffx * z_mult_dat[ip] * (2.0 / dx) * permxp[ip] / viscosity;

              switch (fdir[0])
              {
                case -1:
                  {
                    op = wp;
                    prod_val = rpp[ip - 1] * den_d;
                    diff = value - pp[ip];
                    o_temp = coeff
                             * (diff * RPMean(value, pp[ip], 0.0, prod_der)
                                - RPMean(value, pp[ip], prod_val, prod));
                    break;
                  }

                case 1:
                  {
                    op = ep;
                    prod_val = rpp[ip + 1] * den_d;
                    diff = pp[ip] - value;
                    o_temp = -coeff
                             * (diff * RPMean(pp[ip], value, prod_der, 0.0)
                                + RPMean(pp[ip], value, prod, prod_val));
                    break;
                  }
              }          /* End switch on fdir[0] */
            }            /* End if (fdir[0]) */

            else if (fdir[1])
            {
              coeff = dt * ffy * z_mult_dat[ip] * (2.0 / dy) * permyp[ip] / viscosity;

              switch (fdir[1])
              {
                case -1:
                  {
                    op = sop;
                    prod_val = rpp[ip - sy_v] * den_d;
                    diff = value - pp[ip];
                    o_temp = coeff
                             * (diff * RPMean(value, pp[ip], 0.0, prod_der)
                                - RPMean(value, pp[ip], prod_val, prod));
                    break;
                  }

                case 1:
                  {
                    op = np;
                    prod_val = rpp[ip + sy_v] * den_d;
                    diff = pp[ip] - value;
                    o_temp = -coeff
                             * (diff * RPMean(pp[ip], value, prod_der, 0.0)
                                + RPMean(pp[ip], value, prod, prod_val));
                    break;
                  }
              }          /* End switch on fdir[1] */
            }            /* End if (fdir[1]) */

            else if (fdir[2])
            {
              coeff = dt * ffz * (2.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]))) * permzp[ip] / viscosity;

              switch (fdir[2])
              {
                case -1:
                  {
                    op = lp;
                    prod_val = rpp[ip - sz_v] * den_d;

                    lower_cond = (value) - 0.5 * dz * z_mult_dat[ip] * den_d * gravity;
                    upper_cond = (pp[ip]) + 0.5 * dz * z_mult_dat[ip] * dp[ip] * gravity;
                    diff = lower_cond - upper_cond;

//                    o_temp = coeff
//                             * (diff * RPMean(lower_cond, upper_cond, 0.0, prod_der)
//                                + ((-1.0 - gravity * 0.5 * dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v]) * ddp[ip])
//                                   * RPMean(lower_cond, upper_cond, prod_val, prod)));

                    o_temp = coeff
                             * (diff * RPMean(lower_cond, upper_cond, 0.0, prod_der)
                                + ((-1.0 - gravity * 0.5 * dz * z_mult_dat[ip] * ddp[ip])
                                   * RPMean(lower_cond, upper_cond, prod_val, prod)));

//printf("jacobian lower BC: o_temp=%f prod_der=%f op=%f \n",o_temp, prod_der, op);

                    break;
                  }

                case 1:
                  {
                    op = up;
                    prod_val = rpp[ip + sz_v] * den_d;

                    lower_cond = (pp[ip]) - 0.5 * dz * z_mult_dat[ip] * dp[ip] * gravity;
                    upper_cond = (value) + 0.5 * dz * z_mult_dat[ip] * den_d * gravity;
                    diff = lower_cond - upper_cond;

                    o_temp = -coeff * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)
                                       + ((1.0 - gravity * 0.5 * dz * z_mult_dat[ip] * ddp[ip])
                                          * RPMean(lower_cond, upper_cond, prod, prod_val)));


                    break;
                  }
              }          /* End switch on fdir[2] */
            }         /* End if (fdir[2]) */

            cp[im] += op[im];
            cp[im] -= o_temp;
//            printf("jacobian: cp[im]+=%f cp[im]-=%f \n",op[im], o_temp);

            op[im] = 0.0;
          });

          break;
        }               /* End DirichletBC case */

        case FluxBC:
        {
          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
          {
            im = SubmatrixEltIndex(J_sub, i, j, k);

            if (fdir[0] == -1)
              op = wp;
            if (fdir[0] == 1)
              op = ep;
            if (fdir[1] == -1)
              op = sop;
            if (fdir[1] == 1)
              op = np;
            if (fdir[2] == -1)
              op = lp;
            if (fdir[2] == 1)
              op = up;

            cp[im] += op[im];
            op[im] = 0.0;
          });

          break;
        }         /* End fluxbc case */

        case OverlandBC:     //sk
        {
          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
          {
            im = SubmatrixEltIndex(J_sub, i, j, k);

            //remove contributions to this row corresponding to boundary
            if (fdir[0] == -1)
              op = wp;
            else if (fdir[0] == 1)
              op = ep;
            else if (fdir[1] == -1)
              op = sop;
            else if (fdir[1] == 1)
              op = np;
            else if (fdir[2] == -1)
              op = lp;
            else       // (fdir[2] ==  1)
            {
              op = up;
              /* check if overland flow kicks in */
              if (!ovlnd_flag)
              {
                ip = SubvectorEltIndex(p_sub, i, j, k);
                if ((pp[ip]) > 0.0)
                {
                  ovlnd_flag = 1;
                }
              }
            }

            cp[im] += op[im];
            op[im] = 0.0;       //zero out entry in row of Jacobian
          });

          switch (public_xtra->type)
          {
            case no_nonlinear_jacobian:
            case not_set:
            {
              assert(1);
            }

            case simple:
            {
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
              {
                if (fdir[2] == 1)
                {
                  ip = SubvectorEltIndex(p_sub, i, j, k);
                  io = SubvectorEltIndex(p_sub, i, j, 0);
                  im = SubmatrixEltIndex(J_sub, i, j, k);

                  if ((pp[ip]) > 0.0)
                  {
                    cp[im] += (vol * z_mult_dat[ip]) / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])) * (dt + 1);
                  }
                }
              });
              break;
            }

            case overland_flow:
            {
              /* Get overland flow contributions - DOK*/
              // SGS can we skip this invocation if !overland_flow?
              //PFModuleInvokeType(OverlandFlowEvalInvoke, overlandflow_module,
              //                (grid, is, bc_struct, ipatch, problem_data, pressure,
              //                 ke_der, kw_der, kn_der, ks_der, NULL, NULL, CALCDER));

              if (overlandspinup == 1)
              {
                /* add flux loss equal to excess head  that overwrites the prior overland flux */
                BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
                {
                  if (fdir[2] == 1)
                  {
                    ip = SubvectorEltIndex(p_sub, i, j, k);
                    io = SubvectorEltIndex(p_sub, i, j, 0);
                    im = SubmatrixEltIndex(J_sub, i, j, k);
                    vol = dx * dy * dz;

                    if ((pp[ip]) >= 0.0)
                    {
                      cp[im] += (vol / dz) * dt * (1.0 + 0.0);                     //LEC
//                      printf("Jac SU: CP=%f im=%d  \n", cp[im], im);
                    }
                    else
                    {
                      cp[im] += 0.0;
                    }
                  }
                });
              }
              else
              {
                /* Get overland flow contributions for using kinematic or diffusive - LEC */
                if (diffusive == 0)
                {
                  /*
                  PFModuleInvokeType(OverlandFlowEvalInvoke, overlandflow_module,
                                     (grid, is, bc_struct, ipatch, problem_data, pressure, old_pressure,
                                     ke_der, kw_der, kn_der, ks_der, NULL, NULL, CALCDER));
                  */

                  int sg = is;
                  double *ke_v = ke_der;
                  double *kw_v = kw_der;
                  double *kn_v = kn_der;
                  double *ks_v = ks_der;
                  double *qx_v = NULL;
                  double *qy_v = NULL;

                  int fcn = CALCDER;

                  Vector      *slope_x = ProblemDataTSlopeX(problem_data);
                  Vector      *slope_y = ProblemDataTSlopeY(problem_data);
                  Vector      *mannings = ProblemDataMannings(problem_data);
                  Vector      *top = ProblemDataIndexOfDomainTop(problem_data);

                  Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub;

                  double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp;

                  double xdir, ydir;
                  double q_lo, q_mid, q_hi;
                  double q_v[3];

                  int ival, sy_v, step;
                  int            *fdir;

                  int i, ii, j, k, ip, io, itop;
                  int i1, j1, k1;

                  p_sub = VectorSubvector(pressure, sg);

                  sx_sub = VectorSubvector(slope_x, sg);
                  sy_sub = VectorSubvector(slope_y, sg);
                  mann_sub = VectorSubvector(mannings, sg);
                  top_sub = VectorSubvector(top, sg);

                  pp = SubvectorData(p_sub);

                  sx_dat = SubvectorData(sx_sub);
                  sy_dat = SubvectorData(sy_sub);
                  mann_dat = SubvectorData(mann_sub);
                  top_dat = SubvectorData(top_sub);

                  sy_v = SubvectorNX(top_sub);

                  if (fcn == CALCFCN)
                  {
                    if (qx_v == NULL || qy_v == NULL)  /* do not return velocity fluxes */
                    {
                      BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
                      {
                        if (fdir[2] == 1)
                        {
                          io = SubvectorEltIndex(sx_sub, i, j, 0);
                          itop = SubvectorEltIndex(top_sub, i, j, 0);

                          /* compute east and west faces */
                          /* First initialize velocities, q_v, for inactive region */
                          q_v[0] = 0.0;
                          q_v[1] = 0.0;
                          q_v[2] = 0.0;

                          for (ii = -1; ii < 2; ii++)
                          {
                            k1 = (int)top_dat[itop + ii];
                            if (k1 >= 0)
                            {
                              ip = SubvectorEltIndex(p_sub, (i + ii), j, k1);

                              if (sx_dat[io + ii] > 0.0)
                                xdir = -1.0;
                              else if (sx_dat[io + ii] < 0.0)
                                xdir = 1.0;
                              else
                                xdir = 0.0;

                              q_v[ii + 1] = xdir * (RPowerR(fabs(sx_dat[io + ii]), 0.5) / mann_dat[io + ii]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0));
                            }
                          }

                          /* compute kw and ke - NOTE: io is for current cell */
                          kw_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);
                          ke_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);

                          /* compute north and south faces */
                          /* First initialize velocities, q_v, for inactive region */
                          q_v[0] = 0.0;
                          q_v[1] = 0.0;
                          q_v[2] = 0.0;

                          for (ii = -1; ii < 2; ii++)
                          {
                            step = ii * sy_v;
                            k1 = (int)top_dat[itop + step];
                            if (k1 >= 0)
                            {
                              ip = SubvectorEltIndex(p_sub, i, (j + ii), k1);

                              if (sy_dat[io + step] > 0.0)
                                ydir = -1.0;
                              else if (sy_dat[io + step] < 0.0)
                                ydir = 1.0;
                              else
                                ydir = 0.0;

                              q_v[ii + 1] = ydir * (RPowerR(fabs(sy_dat[io + step]), 0.5) / mann_dat[io + step]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0));
                            }
                          }

                          /* compute ks and kn - NOTE: io is for current cell */
                          ks_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);
                          kn_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);
                        }
                      });
                    }
                    else   /* return velocity fluxes */
                    {
                      BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
                      {
                        if (fdir[2] == 1)
                        {
                          io = SubvectorEltIndex(sx_sub, i, j, 0);
                          itop = SubvectorEltIndex(top_sub, i, j, 0);

                          /* compute east and west faces */
                          /* First initialize velocities, q_v, for inactive region */
                          q_v[0] = 0.0;
                          q_v[1] = 0.0;
                          q_v[2] = 0.0;

                          for (ii = -1; ii < 2; ii++)
                          {
                            k1 = (int)top_dat[itop + ii];
                            if (k1 >= 0)
                            {
                              ip = SubvectorEltIndex(p_sub, (i + ii), j, k1);

                              if (sx_dat[io + ii] > 0.0)
                                xdir = -1.0;
                              else if (sx_dat[io + ii] < 0.0)
                                xdir = 1.0;
                              else
                                xdir = 0.0;

                              q_v[ii + 1] = xdir * (RPowerR(fabs(sx_dat[io + ii]), 0.5) / mann_dat[io + ii]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0));
                            }
                          }
                          qx_v[io] = q_v[1];
                          /* compute kw and ke - NOTE: io is for current cell */
                          kw_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);
                          ke_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);

                          /* compute north and south faces */
                          /* First initialize velocities, q_v, for inactive region */
                          q_v[0] = 0.0;
                          q_v[1] = 0.0;
                          q_v[2] = 0.0;

                          for (ii = -1; ii < 2; ii++)
                          {
                            step = ii * sy_v;
                            k1 = (int)top_dat[itop + step];
                            if (k1 >= 0)
                            {
                              ip = SubvectorEltIndex(p_sub, i, (j + ii), k1);

                              if (sy_dat[io + step] > 0.0)
                                ydir = -1.0;
                              else if (sy_dat[io + step] < 0.0)
                                ydir = 1.0;
                              else
                                ydir = 0.0;

                              q_v[ii + 1] = ydir * (RPowerR(fabs(sy_dat[io + step]), 0.5) / mann_dat[io + step]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0));
                            }
                          }
                          qy_v[io] = q_v[1];
                          /* compute ks and kn - NOTE: io is for current cell */
                          ks_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);
                          kn_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);
                        }
                      });
                    }
                  }
                  else  /* fcn == CALCDER: derivs of KE,KW,KN,KS w.r.t. current cell (i,j,k) */
                  {
                    if (qx_v == NULL || qy_v == NULL)  /* Do not return derivs of velocity fluxes */
                    {
                      BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
                      {
                        if (fdir[2] == 1)
                        {
                          /* compute derivs for east and west faces */

                          /* current cell */
                          io = SubvectorEltIndex(sx_sub, i, j, 0);
                          ip = SubvectorEltIndex(p_sub, i, j, k);

                          if (sx_dat[io] > 0.0)
                            xdir = -1.0;
                          else if (sx_dat[io] < 0.0)
                            xdir = 1.0;
                          else
                            xdir = 0.0;

                          q_mid = xdir * (5.0 / 3.0) * (RPowerR(fabs(sx_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0));
                          /* compute derivs of kw and ke - NOTE: io is for current cell */
                          kw_v[io] = -pfmax(-q_mid, 0.0);
                          ke_v[io] = pfmax(q_mid, 0.0);


                          /* compute north and south faces */
                          if (sy_dat[io] > 0.0)
                            ydir = -1.0;
                          else if (sy_dat[io] < 0.0)
                            ydir = 1.0;
                          else
                            ydir = 0.0;

                          q_mid = ydir * (5.0 / 3.0) * (RPowerR(fabs(sy_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0));
                          /* compute derivs of ks and kn - NOTE: io is for current cell */
                          ks_v[io] = -pfmax(-q_mid, 0.0);
                          kn_v[io] = pfmax(q_mid, 0.0);
                        }
                      });
                    }
                    else   /* return derivs of velocity fluxes */
                    {
                      BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
                      {
                        if (fdir[2] == 1)
                        {
                          /* compute derivs for east and west faces */

                          /* current cell */
                          io = SubvectorEltIndex(sx_sub, i, j, 0);
                          ip = SubvectorEltIndex(p_sub, i, j, k);

                          if (sx_dat[io] > 0.0)
                            xdir = -1.0;
                          else if (sx_dat[io] < 0.0)
                            xdir = 1.0;
                          else
                            xdir = 0.0;

                          q_mid = xdir * (5.0 / 3.0) * (RPowerR(fabs(sx_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0));
                          qx_v[io] = q_mid;
                          /* compute derivs of kw and ke - NOTE: io is for current cell */
                          kw_v[io] = -pfmax(-q_mid, 0.0);
                          ke_v[io] = pfmax(q_mid, 0.0);


                          /* compute north and south faces */
                          if (sy_dat[io] > 0.0)
                            ydir = -1.0;
                          else if (sy_dat[io] < 0.0)
                            ydir = 1.0;
                          else
                            ydir = 0.0;

                          q_mid = ydir * (5.0 / 3.0) * (RPowerR(fabs(sy_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0));
                          qy_v[io] = q_mid;
                          /* compute derivs of ks and kn - NOTE: io is for current cell */
                          ks_v[io] = -pfmax(-q_mid, 0.0);
                          kn_v[io] = pfmax(q_mid, 0.0);
                        }
                      });
                    }
                  }

                }
                else
                {
                  /* Test running Diffuisve calc FCN */
                  //double *dummy1, *dummy2, *dummy3, *dummy4;
                  //PFModuleInvokeType(OverlandFlowEvalDiffInvoke, overlandflow_module_diff, (grid, is, bc_struct, ipatch, problem_data, pressure,
                  //                                             ke_der, kw_der, kn_der, ks_der,
                  //       dummy1, dummy2, dummy3, dummy4,
                  //                                                    NULL, NULL, CALCFCN));
                  /*
                    PFModuleInvokeType(OverlandFlowEvalDiffInvoke, overlandflow_module_diff,
                    (grid, is, bc_struct, ipatch, problem_data, pressure, old_pressure,
                    ke_der, kw_der, kn_der, ks_der,
                    kens_der, kwns_der, knns_der, ksns_der, NULL, NULL, CALCDER));
                  */

                  int sg = is;
                  double *ke_v = ke_der;
                  double *kw_v = kw_der;
                  double *kn_v = kn_der;
                  double *ks_v = ks_der;
                  double *ke_vns = kens_der;
                  double *kw_vns = kwns_der;
                  double *kn_vns = knns_der;
                  double *ks_vns = ksns_der;

                  double *qx_v = NULL;
                  double *qy_v = NULL;
                  int fcn = CALCDER;

                  Vector      *slope_x = ProblemDataTSlopeX(problem_data);
                  Vector      *slope_y = ProblemDataTSlopeY(problem_data);
                  Vector      *mannings = ProblemDataMannings(problem_data);
                  Vector      *top = ProblemDataIndexOfDomainTop(problem_data);

                  // printf("overland_eval_diffusive called\n");
                  Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub, *op_sub;

                  Subgrid      *subgrid;

                  double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp, *opp;

                  double xdir, ydir;
                  double q_lo, q_mid, q_hi;
                  double q_v[4], slope_fx_lo, slope_fx_hi, slope_fx_mid;
                  double slope_fy_lo, slope_fy_hi, slope_fy_mid, dx, dy;
                  double coeff, Pmean, P2, P3, Pdel, Pcen;
                  double slope_mean, manning, s1, s2, Sf_mag;
                  double Press_x, Press_y, Sf_x, Sf_y, Sf_xo, Sf_yo;
                  double Pupx, Pupy, Pupox, Pupoy, Pdown, Pdowno;
                  double ov_epsilon;

                  int ival, sy_v, step;
                  int            *fdir;

                  int i, ii, j, k, ip, ip2, ip3, ip4, ip0, io, itop,k1x, k1y, ipp1, ippsy;
                  int i1, j1, k1, k0x, k0y, iojm1, iojp1, ioip1, ioim1;
                  /* @RMM get grid from global (assuming this is comp grid) to pass to CLM */
                  int gnx = BackgroundNX(GlobalsBackground);
                  int gny = BackgroundNY(GlobalsBackground);

                  p_sub = VectorSubvector(pressure, sg);
                  op_sub = VectorSubvector(old_pressure, sg);
                  sx_sub = VectorSubvector(slope_x, sg);
                  sy_sub = VectorSubvector(slope_y, sg);
                  mann_sub = VectorSubvector(mannings, sg);
                  top_sub = VectorSubvector(top, sg);

                  pp = SubvectorData(p_sub);
                  opp = SubvectorData(op_sub);

                  sx_dat = SubvectorData(sx_sub);
                  sy_dat = SubvectorData(sy_sub);
                  mann_dat = SubvectorData(mann_sub);
                  top_dat = SubvectorData(top_sub);

                  subgrid = GridSubgrid(grid, sg);
                  dx = SubgridDX(subgrid);
                  dy = SubgridDY(subgrid);

                  sy_v = SubvectorNX(top_sub);

                  ov_epsilon = GetDoubleDefault("Solver.OverlandDiffusive.Epsilon", 1.0e-5);
                  if (fcn == CALCFCN)
                  {
                    BCStructPatchLoopOvrlnd(i, j, k, fdir, ival, bc_struct, ipatch, sg,
                    {
                      if (fdir[2] == 1)
                      {
                        io = SubvectorEltIndex(sx_sub, i, j, 0);
                        itop = SubvectorEltIndex(top_sub, i, j, 0);

                        k1 = (int)top_dat[itop];
                        k0x = (int)top_dat[itop - 1];
                        k0y = (int)top_dat[itop - sy_v];
                        k1x = (int)top_dat[itop + 1];
                        k1y = (int)top_dat[itop + sy_v];

                        if (k1 >= 0)
                        {
                          ip = SubvectorEltIndex(p_sub, i, j, k1);
                          ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);
                          ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);
                          Pupx = pfmax(pp[ipp1], 0.0);
                          Pupy = pfmax(pp[ippsy], 0.0);
                          Pupox = pfmax(opp[ipp1], 0.0);
                          Pupoy = pfmax(opp[ippsy], 0.0);
                          Pdown = pfmax(pp[ip], 0.0);
                          Pdowno = pfmax(opp[ip], 0.0);

                          Sf_x = sx_dat[io] + (Pupx - Pdown) / dx;
                          Sf_y = sy_dat[io] + (Pupy - Pdown) / dy;

                          Sf_xo = sx_dat[io] + (Pupox - Pdowno) / dx;
                          Sf_yo = sy_dat[io] + (Pupoy - Pdowno) / dy;

                          Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5);
                          if (Sf_mag < ov_epsilon)
                            Sf_mag = ov_epsilon;

                          Press_x = RPMean(-Sf_x, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ipp1]), 0.0));
                          Press_y = RPMean(-Sf_y, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ippsy]), 0.0));

                          qx_v[io] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                          qy_v[io] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));
                        }

                        //fix for lower x boundary
                        if (k0x < 0.0)
                        {
                          Press_x = pfmax((pp[ip]), 0.0);
                          Sf_x = sx_dat[io] + (Press_x - 0.0) / dx;

                          Pupox = pfmax(opp[ip], 0.0);
                          Sf_xo = sx_dat[io] + (Pupox - 0.0) / dx;

                          double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //+ov_epsilon;
                          if (Sf_mag < ov_epsilon)
                            Sf_mag = ov_epsilon;
                          if (Sf_x > 0.0)
                          {
                            qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                          }
                        }

                        //fix for lower y boundary
                        if (k0y < 0.0)
                        {
                          Press_y = pfmax((pp[ip]), 0.0);
                          Sf_y = sy_dat[io] + (Press_y - 0.0) / dx;

                          Pupoy = pfmax(opp[ip], 0.0);
                          Sf_yo = sy_dat[io] + (Pupoy - 0.0) / dx;

                          double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //Note that the sf_xo was already corrected above
                          if (Sf_mag < ov_epsilon)
                            Sf_mag = ov_epsilon;

                          if (Sf_y > 0.0)
                          {
                            qy_v[io - sy_v] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));
                          }

                          // Recalculating the x flow in the case with both the lower and left boundaries
                          // This is exactly the same as the q_x in the left boundary conditional above but
                          // recalculating qx_v here again becuase the sf_mag will be adjusted with the new sf_yo above
                          if (k0x < 0.0)
                          {
                            if (Sf_x > 0.0)
                            {
                              qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                            }
                          }
                        }
                      }
                    });

                    BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
                    {
                      if (fdir[2] == 1)
                      {
                        io = SubvectorEltIndex(sx_sub, i, j, 0);
                        ke_v[io] = qx_v[io];
                        kw_v[io] = qx_v[io - 1];
                        kn_v[io] = qy_v[io];
                        ks_v[io] = qy_v[io - sy_v];
                        //printf("i=%d j=%d k=%d ke_v=%d kw_v=%d kn_v=%d ks_v=%f\n",i,j,k,ke_v[io],kw_v[io],kn_v[io],ks_v[io]);
                      }
                    });
                  }
                  else          //fcn = CALCDER calculates the derivs of KE KW KN KS wrt to current cell (i,j,k)
                  {
                    BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
                    {
                      if (fdir[2] == 1)
                      {
                        io = SubvectorEltIndex(sx_sub, i, j, 0);
                        itop = SubvectorEltIndex(top_sub, i, j, 0);

                        k1 = (int)top_dat[itop];
                        k0x = (int)top_dat[itop - 1];
                        k0y = (int)top_dat[itop - sy_v];
                        k1x = (int)top_dat[itop + 1];
                        k1y = (int)top_dat[itop + sy_v];

                        if (k1 >= 0)
                        {
                          ip = SubvectorEltIndex(p_sub, i, j, k1);
                          ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);
                          ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);
                          Pupx = pfmax(pp[ipp1], 0.0);
                          Pupy = pfmax(pp[ippsy], 0.0);
                          Pupox = pfmax(opp[ipp1], 0.0);
                          Pupoy = pfmax(opp[ippsy], 0.0);
                          Pdown = pfmax(pp[ip], 0.0);
                          Pdowno = pfmax(opp[ip], 0.0);

                          Sf_x = sx_dat[io] + (Pupx - Pdown) / dx;
                          Sf_y = sy_dat[io] + (Pupy - Pdown) / dy;

                          Sf_xo = sx_dat[io] + (Pupox - Pdowno) / dx;
                          Sf_yo = sy_dat[io] + (Pupoy - Pdowno) / dy;

                          Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //+ov_epsilon;
                          if (Sf_mag < ov_epsilon)
                            Sf_mag = ov_epsilon;

                          if (Sf_x < 0)
                          {
                            ke_v[io] = (5.0 / 3.0) * (-sx_dat[io] - (Pupx / dx)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pdown, (2.0 / 3.0)) +
                                       (8.0 / 3.0) * RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                            kw_v[io + 1] = -RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                            ke_vns[io] = kw_v[io + 1];
                            kw_vns[io + 1] = ke_v[io];
                          }

                          if (Sf_x >= 0)
                          {
                            ke_v[io] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                            kw_v[io + 1] = (5.0 / 3.0) * (-sx_dat[io] + (Pdown / dx)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) -
                                           (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                            ke_vns[io] = kw_v[io + 1];
                            kw_vns[io + 1] = ke_v[io];
                          }

                          if (Sf_y < 0)
                          {
                            kn_v[io] = (5.0 / 3.0) * (-sy_dat[io] - (Pupy / dy)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pdown, (2.0 / 3.0)) +
                                       (8.0 / 3.0) * RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                            ks_v[io + sy_v] = -RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                            kn_vns[io] = ks_v[io + sy_v];
                            ks_vns[io + sy_v] = kn_v[io];
                          }

                          if (Sf_y >= 0)
                          {
                            kn_v[io] = RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                            ks_v[io + sy_v] = (5.0 / 3.0) * (-sy_dat[io] + (Pdown / dy)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupy, (2.0 / 3.0)) -
                                              (8.0 / 3.0) * RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                            kn_vns[io] = ks_v[io + sy_v];
                            ks_vns[io + sy_v] = kn_v[io];
                          }
                        }

                        //fix for lower x boundary
                        if (k0x < 0.0)
                        {
                          Pupx = pfmax((pp[ip]), 0.0);
                          Sf_x = sx_dat[io] + (Pupx - 0.0) / dx;

                          Pupox = pfmax(opp[ip], 0.0);
                          Sf_xo = sx_dat[io] + (Pupox - 0.0) / dx;

                          double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5);
                          if (Sf_mag < ov_epsilon)
                            Sf_mag = ov_epsilon;

                          if (Sf_x < 0)
                          {
                            ke_v[io - 1] = 0.0;
                            kw_v[io] = 0.0;
                            kw_vns[io] = 0.0;
                            ke_vns[io - 1] = 0.0;
                          }

                          if (Sf_x >= 0)
                          {
                            ke_v[io - 1] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                            kw_v[io] = (5.0 / 3.0) * (-sx_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) -
                                       (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                            ke_vns[io - 1] = kw_v[io];
                            kw_vns[io] = ke_v[io - 1];
                          }
                        }

                        //fix for lower y boundary
                        if (k0y < 0.0)
                        {
                          Pupy = pfmax((pp[ip]), 0.0);
                          Sf_y = sy_dat[io] + (Pupy - 0.0) / dy;

                          Pupoy = pfmax(opp[ip], 0.0);
                          Sf_yo = sy_dat[io] + (Pupoy - 0.0) / dy;

                          double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //Note that the sf_xo was already corrected above
                          if (Sf_mag < ov_epsilon)
                            Sf_mag = ov_epsilon;

                          if (Sf_y < 0)
                          {
                            kn_v[io - sy_v] = 0.0;
                            ks_v[io] = 0.0;
                            ks_vns[io] = 0.0;
                            kn_vns[io - sy_v] = 0.0;
                          }

                          if (Sf_y >= 0)
                          {
                            kn_vns[io - sy_v] = RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);
                            ks_v[io] = (5.0 / 3.0) * (-sy_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupy, (2.0 / 3.0)) -
                                       (8.0 / 3.0) * RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);
                            kn_vns[io - sy_v] = ks_v[io];
                            ks_vns[io] = kn_v[io - sy_v];
                          }

                          // Recalculating the x flow in the case with both the lower and left boundaries
                          // This is exactly the same as the q_x in the left boundary conditional above but
                          // recalculating qx_v here again becuase the sf_mag will be adjusted with the new sf_yo above
                          if (k0x < 0.0)
                          {
                            if (Sf_x < 0)
                            {
                              kn_v[io - sy_v] = 0.0;
                              ks_v[io] = 0.0;
                              ks_vns[io] = 0.0;
                              kn_vns[io - sy_v] = 0.0;
                            }

                            if (Sf_x >= 0)
                            {
                              ke_v[io - 1] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                              kw_v[io] = (5.0 / 3.0) * (-sx_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) -
                                         (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                              ke_vns[io - 1] = kw_v[io];
                              kw_vns[io] = ke_v[io - 1];
                            }
                          }
                        }
                      }
                    });
                    //}
                  }
                }
              }


              break;
            }
          }
          break;
        }         /* End overland flow case */

        case SeepageFaceBC:
        {
          /* add flux loss equal to excess head  that overwrites the prior overland flux */
          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
          {
            if (fdir[2] == 1)
            {
              ip = SubvectorEltIndex(p_sub, i, j, k);
              io = SubvectorEltIndex(p_sub, i, j, 0);
              im = SubmatrixEltIndex(J_sub, i, j, k);
              vol = dx * dy * dz;

              if ((pp[ip]) >= 0.0)
              {
                cp[im] += (vol / dz) * dt * (1.0 + 0.0);                       //@RMM
//                  printf("Jac SF: CP=%f im=%d  \n", cp[im], im);
              }
              else
              {
                cp[im] += 0.0;
              }
            }
          });
          break;
        }  // end case seepage face

        /*  OverlandBC for KWE with upwinding, call module */
        case OverlandKinematicBC:
        {
          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
          {
            im = SubmatrixEltIndex(J_sub, i, j, k);
            //remove contributions to this row corresponding to boundary
            if (fdir[0] == -1)
              op = wp;
            else if (fdir[0] == 1)
              op = ep;
            else if (fdir[1] == -1)
              op = sop;
            else if (fdir[1] == 1)
              op = np;
            else if (fdir[2] == -1)
              op = lp;
            else       // (fdir[2] ==  1)
            {
              op = up;
              /* check if overland flow kicks in */
              if (!ovlnd_flag)
              {
                ip = SubvectorEltIndex(p_sub, i, j, k);
                if ((pp[ip]) > 0.0)
                {
                  ovlnd_flag = 1;
                }
              }
            }

            cp[im] += op[im];
            op[im] = 0.0;       //zero out entry in row of Jacobian
          });
          /*
          PFModuleInvokeType(OverlandFlowEvalKinInvoke, overlandflow_module_kin,
                             (grid, is, bc_struct, ipatch, problem_data, pressure,
                              ke_der, kw_der, kn_der, ks_der,
                              NULL, NULL, NULL, NULL, NULL, NULL, CALCDER));
          */

          int sg = is;
          double *ke_v = ke_der;
          double *kw_v = kw_der;
          double *kn_v = kn_der;
          double *ks_v = ks_der;
          double *ke_vns = NULL;
          double *kw_vns = NULL;
          double *kn_vns = NULL;
          double *ks_vns = NULL;

          double *qx_v = NULL;
          double *qy_v = NULL;
          int fcn = CALCDER;

          Vector      *slope_x = ProblemDataTSlopeX(problem_data);
          Vector      *slope_y = ProblemDataTSlopeY(problem_data);
          Vector      *mannings = ProblemDataMannings(problem_data);
          Vector      *top = ProblemDataIndexOfDomainTop(problem_data);

          Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub;

          Subgrid      *subgrid;

          double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp;

          double xdir, ydir;
          double q_lo, q_mid, q_hi, qx_temp, qy_temp;
          double q_v[4], slope_fx_lo, slope_fx_hi, slope_fx_mid;
          double slope_fy_lo, slope_fy_hi, slope_fy_mid, dx, dy;
          double coeff, Pmean, P2, P3, Pdel, Pcen;
          double slope_mean, manning, s1, s2, Sf_mag;
          double Press_x, Press_y, Sf_x, Sf_y, Sf_xo, Sf_yo;
          double ov_epsilon;

          int ival, sy_v, step;
          int            *fdir;

          int i, ii, j, k, ip, ip2, ip3, ip4, ip0, io, itop, k1x, k1y, ipp1, ippsy;
          int i1, j1, k1, k0x, k0y, iojm1, iojp1, ioip1, ioim1;
          /* @RMM get grid from global (assuming this is comp grid) to pass to CLM */
          int gnx = BackgroundNX(GlobalsBackground);
          int gny = BackgroundNY(GlobalsBackground);

          p_sub = VectorSubvector(pressure, sg);

          sx_sub = VectorSubvector(slope_x, sg);
          sy_sub = VectorSubvector(slope_y, sg);
          mann_sub = VectorSubvector(mannings, sg);
          top_sub = VectorSubvector(top, sg);

          pp = SubvectorData(p_sub);

          sx_dat = SubvectorData(sx_sub);
          sy_dat = SubvectorData(sy_sub);
          mann_dat = SubvectorData(mann_sub);
          top_dat = SubvectorData(top_sub);

          subgrid = GridSubgrid(grid, sg);
          dx = SubgridDX(subgrid);
          dy = SubgridDY(subgrid);

          sy_v = SubvectorNX(top_sub);

          //ov_epsilon= 1.0e-5;
          ov_epsilon = GetDoubleDefault("Solver.OverlandKinematic.Epsilon", 1.0e-5);
          if (fcn == CALCFCN)
          {
            BCStructPatchLoopOvrlnd(i, j, k, fdir, ival, bc_struct, ipatch, sg,
            {
              if (fdir[2] == 1)
              {
                io = SubvectorEltIndex(sx_sub, i, j, 0);
                itop = SubvectorEltIndex(top_sub, i, j, 0);

                k1 = (int)top_dat[itop];
                k0x = (int)top_dat[itop - 1];
                k0y = (int)top_dat[itop - sy_v];
                k1x = pfmax((int)top_dat[itop + 1],0);
                k1y = pfmax((int)top_dat[itop + sy_v],0);

                if (k1 >= 0)
                {
                  ip = SubvectorEltIndex(p_sub, i, j, k1);
                  Sf_x = sx_dat[io];
                  Sf_y = sy_dat[io];
                  ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);
                  ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);

                  Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;

                  Press_x = RPMean(-Sf_x, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ipp1]), 0.0));
                  Press_y = RPMean(-Sf_y, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ippsy]), 0.0));

                  qx_v[io] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                  qy_v[io] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));
                }

                //fix for lower x boundary
                if (k0x < 0.0)
                {
                  if (k1 >= 0.0)
                  {
                    Sf_x = sx_dat[io];
                    Sf_y = sy_dat[io];

                    double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);
                    if (Sf_mag < ov_epsilon)
                      Sf_mag = ov_epsilon;

                    if (Sf_x > 0.0)
                    {
                      ip = SubvectorEltIndex(p_sub, i, j, k1);
                      Press_x = pfmax((pp[ip]), 0.0);
                      qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                    }
                  }
                }

                //fix for lower y boundary
                if (k0y < 0.0)
                {
                  if (k1 >= 0.0)
                  {
                    Sf_x = sx_dat[io];
                    Sf_y = sy_dat[io];

                    double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);
                    if (Sf_mag < ov_epsilon)
                      Sf_mag = ov_epsilon;

                    if (Sf_y > 0.0)
                    {
                      ip = SubvectorEltIndex(p_sub, i, j, k1);
                      Press_y = pfmax((pp[ip]), 0.0);
                      qy_v[io - sy_v] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));
                    }
                  }
                }
              }
            });

            BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
            {
              if (fdir[2] == 1)
              {
                io = SubvectorEltIndex(sx_sub, i, j, 0);
                ke_v[io] = qx_v[io];
                kw_v[io] = qx_v[io - 1];
                kn_v[io] = qy_v[io];
                ks_v[io] = qy_v[io - sy_v];
              }
            });
          }
          else          //fcn = CALCDER calculates the derivs
          {
            BCStructPatchLoopOvrlnd(i, j, k, fdir, ival, bc_struct, ipatch, sg,
            {
              if (fdir[2] == 1)
              {
                io = SubvectorEltIndex(sx_sub, i, j, 0);
                itop = SubvectorEltIndex(top_sub, i, j, 0);

                k1 = (int)top_dat[itop];
                k0x = (int)top_dat[itop - 1];
                k0y = (int)top_dat[itop - sy_v];
                k1x = (int)top_dat[itop + 1];
                k1y = (int)top_dat[itop + sy_v];

                if (k1 >= 0)
                {
                  ip = SubvectorEltIndex(p_sub, i, j, k1);
                  ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);
                  ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);

                  Sf_x = sx_dat[io];
                  Sf_y = sy_dat[io];

                  Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;

                  Press_x = RPMean(-Sf_x, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ipp1]), 0.0));
                  Press_y = RPMean(-Sf_y, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ippsy]), 0.0));

                  qx_temp = -(5.0 / 3.0) * (Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (2.0 / 3.0));
                  qy_temp = -(5.0 / 3.0) * (Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (2.0 / 3.0));

                  ke_v[io] = pfmax(qx_temp, 0);
                  kw_v[io + 1] = -pfmax(-qx_temp, 0);
                  kn_v[io] = pfmax(qy_temp, 0);
                  ks_v[io + sy_v] = -pfmax(-qy_temp, 0);
                }

                //fix for lower x boundary
                if (k0x < 0.0)
                {
                  if (k1 >= 0.0)
                  {
                    Sf_x = sx_dat[io];
                    Sf_y = sy_dat[io];

                    double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);
                    if (Sf_mag < ov_epsilon)
                      Sf_mag = ov_epsilon;

                    if (Sf_x > 0.0)
                    {
                      ip = SubvectorEltIndex(p_sub, i, j, k1);
                      Press_x = pfmax((pp[ip]), 0.0);
                      qx_temp = -(5.0 / 3.0) * (Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (2.0 / 3.0));

                      kw_v[io] = qx_temp;
                      ke_v[io - 1] = qx_temp;
                    }
                  }
                }

                //fix for lower y boundary
                if (k0y < 0.0)
                {
                  if (k1 >= 0.0)
                  {
                    Sf_x = sx_dat[io];
                    Sf_y = sy_dat[io];

                    double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);  //+ov_epsilon;
                    if (Sf_mag < ov_epsilon)
                      Sf_mag = ov_epsilon;

                    if (Sf_y > 0.0)
                    {
                      ip = SubvectorEltIndex(p_sub, i, j, k1);
                      Press_y = pfmax((pp[ip]), 0.0);
                      qy_temp = -(5.0 / 3.0) * (Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (2.0 / 3.0));

                      ks_v[io] = qy_temp;
                      kn_v[io - sy_v] = qy_temp;
                    }
                  }
                }
              }
            });
          }   // else calcder


          break;
        } /* End OverlandKinematicBC */

        /* OverlandDiffusiveBC */
        case OverlandDiffusiveBC:
        {
          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
          {
            im = SubmatrixEltIndex(J_sub, i, j, k);

            //remove contributions to this row corresponding to boundary
            if (fdir[0] == -1)
              op = wp;
            else if (fdir[0] == 1)
              op = ep;
            else if (fdir[1] == -1)
              op = sop;
            else if (fdir[1] == 1)
              op = np;
            else if (fdir[2] == -1)
              op = lp;
            else       // (fdir[2] ==  1)
            {
              op = up;
              /* check if overland flow kicks in */
              if (!ovlnd_flag)
              {
                ip = SubvectorEltIndex(p_sub, i, j, k);
                if ((pp[ip]) > 0.0)
                {
                  ovlnd_flag = 1;
                }
              }
            }

            cp[im] += op[im];
            op[im] = 0.0;       //zero out entry in row of Jacobian
          });

          /*
          PFModuleInvokeType(OverlandFlowEvalDiffInvoke, overlandflow_module_diff,
                             (grid, is, bc_struct, ipatch, problem_data, pressure, old_pressure,
                              ke_der, kw_der, kn_der, ks_der,
                              kens_der, kwns_der, knns_der, ksns_der, NULL, NULL, CALCDER));
          */

          int sg = is;
          double *ke_v = ke_der;
          double *kw_v = kw_der;
          double *kn_v = kn_der;
          double *ks_v = ks_der;
          double *ke_vns = kens_der;
          double *kw_vns = kwns_der;
          double *kn_vns = knns_der;
          double *ks_vns = ksns_der;

          double *qx_v = NULL;
          double *qy_v = NULL;
          int fcn = CALCDER;

          Vector      *slope_x = ProblemDataTSlopeX(problem_data);
          Vector      *slope_y = ProblemDataTSlopeY(problem_data);
          Vector      *mannings = ProblemDataMannings(problem_data);
          Vector      *top = ProblemDataIndexOfDomainTop(problem_data);

          // printf("overland_eval_diffusive called\n");
          Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub, *op_sub;

          Subgrid      *subgrid;

          double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp, *opp;

          double xdir, ydir;
          double q_lo, q_mid, q_hi;
          double q_v[4], slope_fx_lo, slope_fx_hi, slope_fx_mid;
          double slope_fy_lo, slope_fy_hi, slope_fy_mid, dx, dy;
          double coeff, Pmean, P2, P3, Pdel, Pcen;
          double slope_mean, manning, s1, s2, Sf_mag;
          double Press_x, Press_y, Sf_x, Sf_y, Sf_xo, Sf_yo;
          double Pupx, Pupy, Pupox, Pupoy, Pdown, Pdowno;
          double ov_epsilon;

          int ival, sy_v, step;
          int            *fdir;

          int i, ii, j, k, ip, ip2, ip3, ip4, ip0, io, itop,k1x, k1y, ipp1, ippsy;
          int i1, j1, k1, k0x, k0y, iojm1, iojp1, ioip1, ioim1;
          /* @RMM get grid from global (assuming this is comp grid) to pass to CLM */
          int gnx = BackgroundNX(GlobalsBackground);
          int gny = BackgroundNY(GlobalsBackground);

          p_sub = VectorSubvector(pressure, sg);
          op_sub = VectorSubvector(old_pressure, sg);
          sx_sub = VectorSubvector(slope_x, sg);
          sy_sub = VectorSubvector(slope_y, sg);
          mann_sub = VectorSubvector(mannings, sg);
          top_sub = VectorSubvector(top, sg);

          pp = SubvectorData(p_sub);
          opp = SubvectorData(op_sub);

          sx_dat = SubvectorData(sx_sub);
          sy_dat = SubvectorData(sy_sub);
          mann_dat = SubvectorData(mann_sub);
          top_dat = SubvectorData(top_sub);

          subgrid = GridSubgrid(grid, sg);
          dx = SubgridDX(subgrid);
          dy = SubgridDY(subgrid);

          sy_v = SubvectorNX(top_sub);

          ov_epsilon = GetDoubleDefault("Solver.OverlandDiffusive.Epsilon", 1.0e-5);
          if (fcn == CALCFCN)
          {
            BCStructPatchLoopOvrlnd(i, j, k, fdir, ival, bc_struct, ipatch, sg,
            {
              if (fdir[2] == 1)
              {
                io = SubvectorEltIndex(sx_sub, i, j, 0);
                itop = SubvectorEltIndex(top_sub, i, j, 0);

                k1 = (int)top_dat[itop];
                k0x = (int)top_dat[itop - 1];
                k0y = (int)top_dat[itop - sy_v];
                k1x = (int)top_dat[itop + 1];
                k1y = (int)top_dat[itop + sy_v];

                if (k1 >= 0)
                {
                  ip = SubvectorEltIndex(p_sub, i, j, k1);
                  ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);
                  ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);
                  Pupx = pfmax(pp[ipp1], 0.0);
                  Pupy = pfmax(pp[ippsy], 0.0);
                  Pupox = pfmax(opp[ipp1], 0.0);
                  Pupoy = pfmax(opp[ippsy], 0.0);
                  Pdown = pfmax(pp[ip], 0.0);
                  Pdowno = pfmax(opp[ip], 0.0);

                  Sf_x = sx_dat[io] + (Pupx - Pdown) / dx;
                  Sf_y = sy_dat[io] + (Pupy - Pdown) / dy;

                  Sf_xo = sx_dat[io] + (Pupox - Pdowno) / dx;
                  Sf_yo = sy_dat[io] + (Pupoy - Pdowno) / dy;

                  Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5);
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;

                  Press_x = RPMean(-Sf_x, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ipp1]), 0.0));
                  Press_y = RPMean(-Sf_y, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ippsy]), 0.0));

                  qx_v[io] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                  qy_v[io] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));
                }

                //fix for lower x boundary
                if (k0x < 0.0)
                {
                  Press_x = pfmax((pp[ip]), 0.0);
                  Sf_x = sx_dat[io] + (Press_x - 0.0) / dx;

                  Pupox = pfmax(opp[ip], 0.0);
                  Sf_xo = sx_dat[io] + (Pupox - 0.0) / dx;

                  double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //+ov_epsilon;
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;
                  if (Sf_x > 0.0)
                  {
                    qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                  }
                }

                //fix for lower y boundary
                if (k0y < 0.0)
                {
                  Press_y = pfmax((pp[ip]), 0.0);
                  Sf_y = sy_dat[io] + (Press_y - 0.0) / dx;

                  Pupoy = pfmax(opp[ip], 0.0);
                  Sf_yo = sy_dat[io] + (Pupoy - 0.0) / dx;

                  double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //Note that the sf_xo was already corrected above
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;

                  if (Sf_y > 0.0)
                  {
                    qy_v[io - sy_v] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));
                  }

                  // Recalculating the x flow in the case with both the lower and left boundaries
                  // This is exactly the same as the q_x in the left boundary conditional above but
                  // recalculating qx_v here again becuase the sf_mag will be adjusted with the new sf_yo above
                  if (k0x < 0.0)
                  {
                    if (Sf_x > 0.0)
                    {
                      qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                    }
                  }
                }
              }
            });

            BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
            {
              if (fdir[2] == 1)
              {
                io = SubvectorEltIndex(sx_sub, i, j, 0);
                ke_v[io] = qx_v[io];
                kw_v[io] = qx_v[io - 1];
                kn_v[io] = qy_v[io];
                ks_v[io] = qy_v[io - sy_v];
                //printf("i=%d j=%d k=%d ke_v=%d kw_v=%d kn_v=%d ks_v=%f\n",i,j,k,ke_v[io],kw_v[io],kn_v[io],ks_v[io]);
              }
            });
          }
          else          //fcn = CALCDER calculates the derivs of KE KW KN KS wrt to current cell (i,j,k)
          {
            BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
            {
              if (fdir[2] == 1)
              {
                io = SubvectorEltIndex(sx_sub, i, j, 0);
                itop = SubvectorEltIndex(top_sub, i, j, 0);

                k1 = (int)top_dat[itop];
                k0x = (int)top_dat[itop - 1];
                k0y = (int)top_dat[itop - sy_v];
                k1x = (int)top_dat[itop + 1];
                k1y = (int)top_dat[itop + sy_v];

                if (k1 >= 0)
                {
                  ip = SubvectorEltIndex(p_sub, i, j, k1);
                  ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);
                  ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);
                  Pupx = pfmax(pp[ipp1], 0.0);
                  Pupy = pfmax(pp[ippsy], 0.0);
                  Pupox = pfmax(opp[ipp1], 0.0);
                  Pupoy = pfmax(opp[ippsy], 0.0);
                  Pdown = pfmax(pp[ip], 0.0);
                  Pdowno = pfmax(opp[ip], 0.0);

                  Sf_x = sx_dat[io] + (Pupx - Pdown) / dx;
                  Sf_y = sy_dat[io] + (Pupy - Pdown) / dy;

                  Sf_xo = sx_dat[io] + (Pupox - Pdowno) / dx;
                  Sf_yo = sy_dat[io] + (Pupoy - Pdowno) / dy;

                  Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //+ov_epsilon;
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;

                  if (Sf_x < 0)
                  {
                    ke_v[io] = (5.0 / 3.0) * (-sx_dat[io] - (Pupx / dx)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pdown, (2.0 / 3.0)) +
                               (8.0 / 3.0) * RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                    kw_v[io + 1] = -RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                    ke_vns[io] = kw_v[io + 1];
                    kw_vns[io + 1] = ke_v[io];
                  }

                  if (Sf_x >= 0)
                  {
                    ke_v[io] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                    kw_v[io + 1] = (5.0 / 3.0) * (-sx_dat[io] + (Pdown / dx)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) -
                                   (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);

                    ke_vns[io] = kw_v[io + 1];
                    kw_vns[io + 1] = ke_v[io];
                  }

                  if (Sf_y < 0)
                  {
                    kn_v[io] = (5.0 / 3.0) * (-sy_dat[io] - (Pupy / dy)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pdown, (2.0 / 3.0)) +
                               (8.0 / 3.0) * RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                    ks_v[io + sy_v] = -RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                    kn_vns[io] = ks_v[io + sy_v];
                    ks_vns[io + sy_v] = kn_v[io];
                  }

                  if (Sf_y >= 0)
                  {
                    kn_v[io] = RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                    ks_v[io + sy_v] = (5.0 / 3.0) * (-sy_dat[io] + (Pdown / dy)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupy, (2.0 / 3.0)) -
                                      (8.0 / 3.0) * RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);

                    kn_vns[io] = ks_v[io + sy_v];
                    ks_vns[io + sy_v] = kn_v[io];
                  }
                }

                //fix for lower x boundary
                if (k0x < 0.0)
                {
                  Pupx = pfmax((pp[ip]), 0.0);
                  Sf_x = sx_dat[io] + (Pupx - 0.0) / dx;

                  Pupox = pfmax(opp[ip], 0.0);
                  Sf_xo = sx_dat[io] + (Pupox - 0.0) / dx;

                  double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5);
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;

                  if (Sf_x < 0)
                  {
                    ke_v[io - 1] = 0.0;
                    kw_v[io] = 0.0;
                    kw_vns[io] = 0.0;
                    ke_vns[io - 1] = 0.0;
                  }

                  if (Sf_x >= 0)
                  {
                    ke_v[io - 1] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                    kw_v[io] = (5.0 / 3.0) * (-sx_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) -
                               (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                    ke_vns[io - 1] = kw_v[io];
                    kw_vns[io] = ke_v[io - 1];
                  }
                }

                //fix for lower y boundary
                if (k0y < 0.0)
                {
                  Pupy = pfmax((pp[ip]), 0.0);
                  Sf_y = sy_dat[io] + (Pupy - 0.0) / dy;

                  Pupoy = pfmax(opp[ip], 0.0);
                  Sf_yo = sy_dat[io] + (Pupoy - 0.0) / dy;

                  double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); //Note that the sf_xo was already corrected above
                  if (Sf_mag < ov_epsilon)
                    Sf_mag = ov_epsilon;

                  if (Sf_y < 0)
                  {
                    kn_v[io - sy_v] = 0.0;
                    ks_v[io] = 0.0;
                    ks_vns[io] = 0.0;
                    kn_vns[io - sy_v] = 0.0;
                  }

                  if (Sf_y >= 0)
                  {
                    kn_vns[io - sy_v] = RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);
                    ks_v[io] = (5.0 / 3.0) * (-sy_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupy, (2.0 / 3.0)) -
                               (8.0 / 3.0) * RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy);
                    kn_vns[io - sy_v] = ks_v[io];
                    ks_vns[io] = kn_v[io - sy_v];
                  }

                  // Recalculating the x flow in the case with both the lower and left boundaries
                  // This is exactly the same as the q_x in the left boundary conditional above but
                  // recalculating qx_v here again becuase the sf_mag will be adjusted with the new sf_yo above
                  if (k0x < 0.0)
                  {
                    if (Sf_x < 0)
                    {
                      kn_v[io - sy_v] = 0.0;
                      ks_v[io] = 0.0;
                      ks_vns[io] = 0.0;
                      kn_vns[io - sy_v] = 0.0;
                    }

                    if (Sf_x >= 0)
                    {
                      ke_v[io - 1] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                      kw_v[io] = (5.0 / 3.0) * (-sx_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) -
                                 (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx);
                      ke_vns[io - 1] = kw_v[io];
                      kw_vns[io] = ke_v[io - 1];
                    }
                  }
                }
              }
            });
            //}
          }

          break;
        } /* End OverlandDiffusiveBC */
      }        /* End switch BCtype */
    }          /* End ipatch loop */
  }            /* End subgrid loop */

  /*
  // Invokes problem_richards_bc_internal.c:RichardsBCInternal
  PFModuleInvokeType(RichardsBCInternalInvoke, bc_internal, (problem, problem_data, NULL, J, time,
                                                             pressure, CALCDER));
  */

  {
    GetModulePublicXtra(RichardsBCInternal, bc_internal, public_xtra);

    //Problem *    problem = problem;
    //ProblemData *problem_data = problem_data;
    Vector *     f = NULL;
    Matrix *     A = J;
    //double       time = time;
    //Vector *     pressure = pressure;
    int          fcn = CALCDER;

    WellData         *well_data = ProblemDataWellData(problem_data);
    WellDataPhysical *well_data_physical;
    WellDataValue    *well_data_value;

    TimeCycleData    *time_cycle_data;

    int num_conditions = (public_xtra->num_conditions);
    int num_wells, total_num;

    Grid             *grid = VectorGrid(pressure);

    SubgridArray     *internal_bc_subgrids = NULL;

    Subgrid          *subgrid, *subgrid_ind, *new_subgrid;

    Subvector        *p_sub;


    double           *pp;
    double           *internal_bc_conditions = NULL;

    double dx, dy, dz;
    double value;

    int ix, iy, iz;
    int nx, ny, nz;
    int rx, ry, rz;
    int process;

    int i, j, k;
    int grid_index, well, index;
    int cycle_number, interval_number;
    int ip, im;

    /*--------------------------------------------------------------------
     * gridify the internal boundary locations (should be done elsewhere?)
     *--------------------------------------------------------------------*/

    if (num_conditions > 0)
    {
      internal_bc_subgrids = NewSubgridArray();
      internal_bc_conditions = ctalloc(double, num_conditions);

      for (i = 0; i < num_conditions; i++)
      {
        switch ((public_xtra->type[i]))
        {
          case 0:
          {
            GetDummyType(RichardsBCInternal, 0, (public_xtra->data[i]), dummy0);

            ix = IndexSpaceX((dummy0->xlocation), 0);
            iy = IndexSpaceY((dummy0->ylocation), 0);
            iz = IndexSpaceZ((dummy0->zlocation), 0);

            nx = 1;
            ny = 1;
            nz = 1;

            rx = 0;
            ry = 0;
            rz = 0;

            process = amps_Rank(amps_CommWorld);

            new_subgrid = NewSubgrid(ix, iy, iz,
                                     nx, ny, nz,
                                     rx, ry, rz,
                                     process);

            AppendSubgrid(new_subgrid, internal_bc_subgrids);

            internal_bc_conditions[i] = (dummy0->value);

            break;
          }
        }
      }
    }

    /*--------------------------------------------------------------------
     * Put in the internal conditions using the subgrids computed above
     * Put in any pressure wells from the well package
     *--------------------------------------------------------------------*/

    num_wells = WellDataNumPressWells(well_data);
    total_num = num_conditions + num_wells;

    if ((num_conditions > 0) || (num_wells > 0))
    {
      /* Set explicit pressure assignments*/

      for (grid_index = 0; grid_index < GridNumSubgrids(grid); grid_index++)
      {
        subgrid = GridSubgrid(grid, grid_index);

        p_sub = VectorSubvector(pressure, grid_index);
        pp = SubvectorData(p_sub);


        for (index = 0; index < total_num; index++)
        {
          if (index < num_conditions)
          {
            subgrid_ind = SubgridArraySubgrid(internal_bc_subgrids, index);
            value = internal_bc_conditions[index];
          }
          else
          {
            well = index - num_conditions;
            time_cycle_data = WellDataTimeCycleData(well_data);
            well_data_physical = WellDataPressWellPhysical(well_data, well);
            cycle_number = WellDataPhysicalCycleNumber(well_data_physical);
            interval_number =
                              TimeCycleDataComputeIntervalNumber(problem, time,
                                                                 time_cycle_data,
                                                                 cycle_number);
            well_data_value =
                              WellDataPressWellIntervalValue(well_data, well,
                                                             interval_number);
            subgrid_ind = WellDataPhysicalSubgrid(well_data_physical);
            value = WellDataValuePhaseValue(well_data_value, 0);
          }

          ix = SubgridIX(subgrid_ind);
          iy = SubgridIY(subgrid_ind);
          iz = SubgridIZ(subgrid_ind);

          nx = SubgridNX(subgrid_ind);
          ny = SubgridNY(subgrid_ind);
          nz = SubgridNZ(subgrid_ind);

          dx = SubgridDX(subgrid_ind);
          dy = SubgridDY(subgrid_ind);
          dz = SubgridDZ(subgrid_ind);

          if (fcn == CALCFCN)
          {
            Subvector *f_sub = VectorSubvector(f, grid_index);
            double *fp = SubvectorData(f_sub);

            BoxLoopI0(i, j, k,
                      ix, iy, iz, nx, ny, nz,
            {
              /* Need to check if i,j,k is part of this subgrid or not */
              if (((i >= SubgridIX(subgrid)) &&
                   (i < SubgridIX(subgrid) + SubgridNX(subgrid))) &&
                  ((j >= SubgridIY(subgrid)) &&
                   (j < SubgridIY(subgrid) + SubgridNY(subgrid))) &&
                  ((k >= SubgridIZ(subgrid)) &&
                   (k < SubgridIZ(subgrid) + SubgridNZ(subgrid))))
              {
                ip = SubvectorEltIndex(f_sub, i, j, k);
                fp[ip] = pp[ip] - value;
              }
            });
          }
          else if (fcn == CALCDER)
          {
            Submatrix        *A_sub = MatrixSubmatrix(A, grid_index);

            double *cp = SubmatrixStencilData(A_sub, 0);
            double *wp = SubmatrixStencilData(A_sub, 1);
            double *ep = SubmatrixStencilData(A_sub, 2);
            double *sp = SubmatrixStencilData(A_sub, 3);
            double *np = SubmatrixStencilData(A_sub, 4);
            double *lp = SubmatrixStencilData(A_sub, 5);
            double *up = SubmatrixStencilData(A_sub, 6);

            BoxLoopI0(i, j, k,
                      ix, iy, iz, nx, ny, nz,
            {
              /* Need to check if i,j,k is part of this subgrid or not */
              if (((i >= SubgridIX(subgrid)) &&
                   (i < SubgridIX(subgrid) + SubgridNX(subgrid))) &&
                  ((j >= SubgridIY(subgrid)) &&
                   (j < SubgridIY(subgrid) + SubgridNY(subgrid))) &&
                  ((k >= SubgridIZ(subgrid)) &&
                   (k < SubgridIZ(subgrid) + SubgridNZ(subgrid))))
              {
                im = SubmatrixEltIndex(A_sub, i, j, k);
                cp[im] = 1.0;
                wp[im] = 0.0;
                ep[im] = 0.0;
                sp[im] = 0.0;
                np[im] = 0.0;
                lp[im] = 0.0;
                up[im] = 0.0;
              }
            });
          }
        }           /* End loop over conditions */
      }             /* End loop over processor subgrids */


      if (num_conditions > 0)
      {
        FreeSubgridArray(internal_bc_subgrids);
      }
      tfree(internal_bc_conditions);
    }               /* End if have well or internal pressure conditions */

  }


  if (public_xtra->type == overland_flow)
  {
    // SGS always have to do communication here since
    // each processor may/may not be doing overland flow.
    /* Update ghost points for JB before building JC */
    if (MatrixCommPkg(J))
    {
      handle = InitMatrixUpdate(J);
      FinalizeMatrixUpdate(handle);
    }

    /* Pass KW values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KW, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
    /* Pass KE values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KE, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
    /* Pass KS values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KS, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
    /* Pass KN values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KN, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
    /* Pass KWns values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KWns, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
    /* Pass KEns values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KEns, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
    /* Pass KSns values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KSns, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
    /* Pass KNns values to neighbors.  */
    vector_update_handle = InitVectorUpdate(KNns, VectorUpdateAll);
    FinalizeVectorUpdate(vector_update_handle);
  }

  /* Build submatrix JC if overland flow case */
  if (ovlnd_flag && public_xtra->type == overland_flow)
  {
    /* begin loop to build JC */
    ForSubgridI(is, GridSubgrids(grid))
    {
      subgrid = GridSubgrid(grid, is);

      dx = SubgridDX(subgrid);
      dy = SubgridDY(subgrid);
      dz = SubgridDZ(subgrid);

      vol = dx * dy * dz;

      ffx = dy * dz;
      ffy = dx * dz;
      ffz = dx * dy;

      p_sub = VectorSubvector(pressure, is);

      J_sub = MatrixSubmatrix(J, is);
      JC_sub = MatrixSubmatrix(JC, is);

      kw_sub = VectorSubvector(KW, is);
      ke_sub = VectorSubvector(KE, is);
      kn_sub = VectorSubvector(KN, is);
      ks_sub = VectorSubvector(KS, is);
      kwns_sub = VectorSubvector(KWns, is);
      kens_sub = VectorSubvector(KEns, is);
      knns_sub = VectorSubvector(KNns, is);
      ksns_sub = VectorSubvector(KSns, is);

      top_sub = VectorSubvector(top, is);
      sx_sub = VectorSubvector(slope_x, is);

      sy_v = SubvectorNX(sx_sub);
      nx_m = SubmatrixNX(J_sub);
      ny_m = SubmatrixNY(J_sub);
      sy_m = nx_m;
      sz_m = nx_m * ny_m;

      ix = SubgridIX(subgrid);
      iy = SubgridIY(subgrid);
      iz = SubgridIZ(subgrid);

      nx = SubgridNX(subgrid);
      ny = SubgridNY(subgrid);

      pp = SubvectorData(p_sub);
      /* for Bmat */
      cp = SubmatrixStencilData(J_sub, 0);
      wp = SubmatrixStencilData(J_sub, 1);
      ep = SubmatrixStencilData(J_sub, 2);
      sop = SubmatrixStencilData(J_sub, 3);
      np = SubmatrixStencilData(J_sub, 4);
      lp = SubmatrixStencilData(J_sub, 5);
      up = SubmatrixStencilData(J_sub, 6);

      /* for Cmat */
      cp_c = SubmatrixStencilData(JC_sub, 0);
      wp_c = SubmatrixStencilData(JC_sub, 1);
      ep_c = SubmatrixStencilData(JC_sub, 2);
      sop_c = SubmatrixStencilData(JC_sub, 3);
      np_c = SubmatrixStencilData(JC_sub, 4);

      kw_der = SubvectorData(kw_sub);
      ke_der = SubvectorData(ke_sub);
      kn_der = SubvectorData(kn_sub);
      ks_der = SubvectorData(ks_sub);
      kwns_der = SubvectorData(kwns_sub);
      kens_der = SubvectorData(kens_sub);
      knns_der = SubvectorData(knns_sub);
      ksns_der = SubvectorData(ksns_sub);

      top_dat = SubvectorData(top_sub);

      for (ipatch = 0; ipatch < BCStructNumPatches(bc_struct); ipatch++)
      {
        switch (BCStructBCType(bc_struct, ipatch))
        {
          /* Fall through cases for new Overland types */
          case OverlandKinematicBC:
          {
            BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
            {
              if (fdir[2] == 1)
              {
                /* Loop over boundary patches to build JC matrix.
                 */
                io = SubmatrixEltIndex(J_sub, i, j, iz);
                io1 = SubvectorEltIndex(sx_sub, i, j, 0);
                itop = SubvectorEltIndex(top_sub, i, j, 0);

                /* Update JC */
                ip = SubvectorEltIndex(p_sub, i, j, k);
                im = SubmatrixEltIndex(J_sub, i, j, k);

                /* First put contributions from subsurface diagonal onto diagonal of JC */
                cp_c[io] = cp[im];
                cp[im] = 0.0;         // update JB
                /* Now check off-diagonal nodes to see if any surface-surface connections exist */
                /* West */
                k1 = (int)top_dat[itop - 1];

                if (k1 >= 0)
                {
                  if (k1 == k)         /*west node is also surface node */
                  {
                    wp_c[io] += wp[im];
                    wp[im] = 0.0;           // update JB
                  }
                }
                /* East */
                k1 = (int)top_dat[itop + 1];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*east node is also surface node */
                  {
                    ep_c[io] += ep[im];
                    ep[im] = 0.0;           //update JB
                  }
                }
                /* South */
                k1 = (int)top_dat[itop - sy_v];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*south node is also surface node */
                  {
                    sop_c[io] += sop[im];
                    sop[im] = 0.0;           //update JB
                  }
                }
                /* North */
                k1 = (int)top_dat[itop + sy_v];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*north node is also surface node */
                  {
                    np_c[io] += np[im];
                    np[im] = 0.0;           // Update JB
                  }
                }

                /* Now add overland contributions to JC */
                if ((pp[ip]) > 0.0)
                {
                  /*diagonal term */
                  cp_c[io] += (vol / dz) + (vol / ffy) * dt * (ke_der[io1] - kw_der[io1])
                              + (vol / ffx) * dt * (kn_der[io1] - ks_der[io1]);
                }

                /*west term */
                wp_c[io] -= (vol / ffy) * dt * (ke_der[io1 - 1]);

                /*East term */
                ep_c[io] += (vol / ffy) * dt * (kw_der[io1 + 1]);

                /*south term */
                sop_c[io] -= (vol / ffx) * dt * (kn_der[io1 - sy_v]);

                /*north term */
                np_c[io] += (vol / ffx) * dt * (ks_der[io1 + sy_v]);
              }
            });
            break;
          }

          case OverlandDiffusiveBC:
          {
            BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
            {
              if (fdir[2] == 1)
              {
                /* Loop over boundary patches to build JC matrix.
                 */
                io = SubmatrixEltIndex(J_sub, i, j, iz);
                io1 = SubvectorEltIndex(sx_sub, i, j, 0);
                itop = SubvectorEltIndex(top_sub, i, j, 0);

                /* Update JC */
                ip = SubvectorEltIndex(p_sub, i, j, k);
                im = SubmatrixEltIndex(J_sub, i, j, k);

                /* First put contributions from subsurface diagonal onto diagonal of JC */
                cp_c[io] = cp[im];
                cp[im] = 0.0;         // update JB
                /* Now check off-diagonal nodes to see if any surface-surface connections exist */
                /* West */
                k1 = (int)top_dat[itop - 1];

                if (k1 >= 0)
                {
                  if (k1 == k)         /*west node is also surface node */
                  {
                    wp_c[io] += wp[im];
                    wp[im] = 0.0;           // update JB
                  }
                }
                /* East */
                k1 = (int)top_dat[itop + 1];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*east node is also surface node */
                  {
                    ep_c[io] += ep[im];
                    ep[im] = 0.0;           //update JB
                  }
                }
                /* South */
                k1 = (int)top_dat[itop - sy_v];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*south node is also surface node */
                  {
                    sop_c[io] += sop[im];
                    sop[im] = 0.0;           //update JB
                  }
                }
                /* North */
                k1 = (int)top_dat[itop + sy_v];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*north node is also surface node */
                  {
                    np_c[io] += np[im];
                    np[im] = 0.0;           // Update JB
                  }
                }

                /* Now add overland contributions to JC */
                if ((pp[ip]) > 0.0)
                {
                  /*diagonal term */
                  cp_c[io] += (vol / dz) + (vol / ffy) * dt * (ke_der[io1] - kw_der[io1])
                              + (vol / ffx) * dt * (kn_der[io1] - ks_der[io1]);
                }
                /*west term */
                wp_c[io] -= (vol / ffy) * dt * (kwns_der[io1]);

                /*East term */
                ep_c[io] += (vol / ffy) * dt * (kens_der[io1]);

                /*south term */
                sop_c[io] -= (vol / ffx) * dt * (ksns_der[io1]);

                /*north term */
                np_c[io] += (vol / ffx) * dt * (knns_der[io1]);
              }
            });
            break;
          }

          case OverlandBC:
          {
            BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,
            {
              if (fdir[2] == 1)
              {
                /* Loop over boundary patches to build JC matrix.
                 */
                io = SubmatrixEltIndex(J_sub, i, j, iz);
                io1 = SubvectorEltIndex(sx_sub, i, j, 0);
                itop = SubvectorEltIndex(top_sub, i, j, 0);

                /* Update JC */
                ip = SubvectorEltIndex(p_sub, i, j, k);
                im = SubmatrixEltIndex(J_sub, i, j, k);

                /* First put contributions from subsurface diagonal onto diagonal of JC */
                cp_c[io] = cp[im];
                cp[im] = 0.0;         // update JB
                /* Now check off-diagonal nodes to see if any surface-surface connections exist */
                /* West */
                k1 = (int)top_dat[itop - 1];

                if (k1 >= 0)
                {
                  if (k1 == k)         /*west node is also surface node */
                  {
                    wp_c[io] += wp[im];
                    wp[im] = 0.0;           // update JB
                  }
                }
                /* East */
                k1 = (int)top_dat[itop + 1];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*east node is also surface node */
                  {
                    ep_c[io] += ep[im];
                    ep[im] = 0.0;           //update JB
                  }
                }
                /* South */
                k1 = (int)top_dat[itop - sy_v];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*south node is also surface node */
                  {
                    sop_c[io] += sop[im];
                    sop[im] = 0.0;           //update JB
                  }
                }
                /* North */
                k1 = (int)top_dat[itop + sy_v];
                if (k1 >= 0)
                {
                  if (k1 == k)         /*north node is also surface node */
                  {
                    np_c[io] += np[im];
                    np[im] = 0.0;           // Update JB
                  }
                }

                /* Now add overland contributions to JC */
                if ((pp[ip]) > 0.0)
                {
                  /*diagonal term */
                  cp_c[io] += (vol / dz) + (vol / ffy) * dt * (ke_der[io1] - kw_der[io1])
                              + (vol / ffx) * dt * (kn_der[io1] - ks_der[io1]);
                }
                else
                {
                  // Laura's version
                  cp_c[io] += 0.0 + dt * (vol / dz) * (public_xtra->SpinupDampP1 * exp(pfmin(pp[ip], 0.0) * public_xtra->SpinupDampP1) * public_xtra->SpinupDampP2); //NBE
                }

                if (diffusive == 0)
                {
                  /*west term */
                  wp_c[io] -= (vol / ffy) * dt * (ke_der[io1 - 1]);

                  /*East term */
                  ep_c[io] += (vol / ffy) * dt * (kw_der[io1 + 1]);

                  /*south term */
                  sop_c[io] -= (vol / ffx) * dt * (kn_der[io1 - sy_v]);

                  /*north term */
                  np_c[io] += (vol / ffx) * dt * (ks_der[io1 + sy_v]);
                }
                else
                {
                  /*west term */
                  wp_c[io] -= (vol / ffy) * dt * (kwns_der[io1]);

                  /*East term */
                  ep_c[io] += (vol / ffy) * dt * (kens_der[io1]);

                  /*south term */
                  sop_c[io] -= (vol / ffx) * dt * (ksns_der[io1]);

                  /*north term */
                  np_c[io] += (vol / ffx) * dt * (knns_der[io1]);
                }
              }
            });
            break;
          }
        }         /* End switch BCtype */
      }           /* End ipatch loop */
    }             /* End subgrid loop */
  }



  /* Set pressures outside domain to zero.
   * Recall: equation to solve is f = 0, so components of f outside
   * domain are set to the respective pressure value.
   *
   * Should change this to set pressures to scaling value.
   * CSW: Should I set this to pressure * vol * dt ??? */

  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    J_sub = MatrixSubmatrix(J, is);

    /* RDF: assumes resolutions are the same in all 3 directions */
    r = SubgridRX(subgrid);

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);
    iz = SubgridIZ(subgrid);

    nx = SubgridNX(subgrid);
    ny = SubgridNY(subgrid);
    nz = SubgridNZ(subgrid);

    cp = SubmatrixStencilData(J_sub, 0);
    wp = SubmatrixStencilData(J_sub, 1);
    ep = SubmatrixStencilData(J_sub, 2);
    sop = SubmatrixStencilData(J_sub, 3);
    np = SubmatrixStencilData(J_sub, 4);
    lp = SubmatrixStencilData(J_sub, 5);
    up = SubmatrixStencilData(J_sub, 6);

    /* for Cmat */
    JC_sub = MatrixSubmatrix(JC, is);
    cp_c = SubmatrixStencilData(JC_sub, 0);
    wp_c = SubmatrixStencilData(JC_sub, 1);
    ep_c = SubmatrixStencilData(JC_sub, 2);
    sop_c = SubmatrixStencilData(JC_sub, 3);
    np_c = SubmatrixStencilData(JC_sub, 4);

    GrGeomOutLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      im = SubmatrixEltIndex(J_sub, i, j, k);
      cp[im] = 1.0;
      wp[im] = 0.0;
      ep[im] = 0.0;
      sop[im] = 0.0;
      np[im] = 0.0;
      lp[im] = 0.0;
      up[im] = 0.0;

//#if 0
//                     /* JC matrix */
//                     cp_c[im] = 1.0;
//                     wp_c[im] = 0.0;
//                     ep_c[im] = 0.0;
//                     sop_c[im] = 0.0;
//                     np_c[im] = 0.0;
//#endif */
    });
  }


  /*-----------------------------------------------------------------------
   * Update matrix ghost points
   *-----------------------------------------------------------------------*/
  if (public_xtra->type == overland_flow)
  {
    /* Update matrices and setup pointers */
    if (MatrixCommPkg(J))
    {
      handle = InitMatrixUpdate(J);
      FinalizeMatrixUpdate(handle);
    }
    *ptr_to_J = J;

    if (MatrixCommPkg(JC))
    {
      handle = InitMatrixUpdate(JC);
      FinalizeMatrixUpdate(handle);
    }
    *ptr_to_JC = JC;
  }
  else  /* No overland flow */
  {
    *ptr_to_JC = NULL;

    if (MatrixCommPkg(J))
    {
      handle = InitMatrixUpdate(J);
      FinalizeMatrixUpdate(handle);
    }

    *ptr_to_J = J;
  } /* end if ovlnd_flag */

  /*-----------------------------------------------------------------------
   * Free temp vectors
   *-----------------------------------------------------------------------*/

  FreeBCStruct(bc_struct);

  FreeVector(density_der);
  FreeVector(saturation_der);
  FreeVector(KW);
  FreeVector(KE);
  FreeVector(KN);
  FreeVector(KS);
  FreeVector(KWns);
  FreeVector(KEns);
  FreeVector(KNns);
  FreeVector(KSns);

  return;
}


/*--------------------------------------------------------------------------
 * RichardsJacobianEvalInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule    *RichardsJacobianEvalInitInstanceXtra(
                                                  Problem *    problem,
                                                  Grid *       grid,
                                                  ProblemData *problem_data,
                                                  double *     temp_data,
                                                  int          symmetric_jac)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra;

  Stencil       *stencil, *stencil_C;

  (void)problem_data;

  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (grid != NULL)
  {
    /* free old data */
    if ((instance_xtra->grid) != NULL)
    {
      FreeMatrix(instance_xtra->J);
      FreeMatrix(instance_xtra->JC);      /* DOK */
    }

    /* set new data */
    (instance_xtra->grid) = grid;

    /* set up jacobian matrix */
    stencil = NewStencil(jacobian_stencil_shape, 7);
    stencil_C = NewStencil(jacobian_stencil_shape, 7);  //DOK

    if (symmetric_jac)
    {
      (instance_xtra->J) = NewMatrixType(grid, NULL, stencil, ON, stencil,
                                         matrix_cell_centered);
      (instance_xtra->JC) = NewMatrixType(grid, NULL, stencil_C, ON, stencil_C,
                                          matrix_cell_centered);
    }
    else
    {
      (instance_xtra->J) = NewMatrixType(grid, NULL, stencil, OFF, stencil,
                                         matrix_cell_centered);
      (instance_xtra->JC) = NewMatrixType(grid, NULL, stencil_C, OFF, stencil_C,
                                          matrix_cell_centered);
    }
  }

  if (temp_data != NULL)
  {
    (instance_xtra->temp_data) = temp_data;
  }

  if (problem != NULL)
  {
    (instance_xtra->problem) = problem;
  }

  if (PFModuleInstanceXtra(this_module) == NULL)
  {
    (instance_xtra->density_module) =
      PFModuleNewInstance(ProblemPhaseDensity(problem), ());
    (instance_xtra->bc_pressure) =
      PFModuleNewInstanceType(BCPressureInitInstanceXtraInvoke, ProblemBCPressure(problem), (problem));
    (instance_xtra->saturation_module) =
      PFModuleNewInstanceType(SaturationInitInstanceXtraInvoke, ProblemSaturation(problem), (NULL, NULL));
    (instance_xtra->rel_perm_module) =
      PFModuleNewInstanceType(PhaseRelPermInitInstanceXtraInvoke, ProblemPhaseRelPerm(problem), (NULL, NULL));
    (instance_xtra->bc_internal) =
      PFModuleNewInstance(ProblemBCInternal(problem), ());
    (instance_xtra->overlandflow_module) =
      PFModuleNewInstance(ProblemOverlandFlowEval(problem), ());     //DOK
    (instance_xtra->overlandflow_module_diff) =
      PFModuleNewInstance(ProblemOverlandFlowEvalDiff(problem), ());   //RMM-LEC
    (instance_xtra->overlandflow_module_kin)
      = PFModuleNewInstance(ProblemOverlandFlowEvalKin(problem), ());
  }
  else
  {
    PFModuleReNewInstance((instance_xtra->density_module), ());
    PFModuleReNewInstanceType(BCPressureInitInstanceXtraInvoke, (instance_xtra->bc_pressure), (problem));
    PFModuleReNewInstanceType(SaturationInitInstanceXtraInvoke, (instance_xtra->saturation_module),
                              (NULL, NULL));
    PFModuleReNewInstanceType(PhaseRelPermInitInstanceXtraInvoke, (instance_xtra->rel_perm_module),
                              (NULL, NULL));
    PFModuleReNewInstance((instance_xtra->bc_internal), ());
    PFModuleReNewInstance((instance_xtra->overlandflow_module), ());     //DOK
    PFModuleReNewInstance((instance_xtra->overlandflow_module_diff), ());      //RMM-LEC
    PFModuleReNewInstance((instance_xtra->overlandflow_module_kin), ());
  }


  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * RichardsJacobianEvalFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  RichardsJacobianEvalFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    PFModuleFreeInstance(instance_xtra->density_module);
    PFModuleFreeInstance(instance_xtra->bc_pressure);
    PFModuleFreeInstance(instance_xtra->saturation_module);
    PFModuleFreeInstance(instance_xtra->rel_perm_module);
    PFModuleFreeInstance(instance_xtra->bc_internal);
    PFModuleFreeInstance(instance_xtra->overlandflow_module);     //DOK
    PFModuleFreeInstance(instance_xtra->overlandflow_module_diff);       //RMM-LEC
    PFModuleFreeInstance(instance_xtra->overlandflow_module_kin);

    FreeMatrix(instance_xtra->J);

    FreeMatrix(instance_xtra->JC);     /* DOK */

    tfree(instance_xtra);
  }
}


/*--------------------------------------------------------------------------
 * RichardsJacobianEvalNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule   *RichardsJacobianEvalNewPublicXtra(char *name)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;
  char key[IDB_MAX_KEY_LEN];
  char          *switch_name;
  int switch_value;
  NameArray switch_na;
  NameArray upwind_switch_na;


  (void)name;

  public_xtra = ctalloc(PublicXtra, 1);
/* These parameters dampen the transition/switching into overland flow to speedup
 * the spinup process */
  sprintf(key, "OverlandSpinupDampP1");
  public_xtra->SpinupDampP1 = GetDoubleDefault(key, 0.0);
  sprintf(key, "OverlandSpinupDampP2");
  public_xtra->SpinupDampP2 = GetDoubleDefault(key, 0.0);    // NBE

  ///* parameters for upwinding formulation for TFG */
  upwind_switch_na = NA_NewNameArray("Original UpwindSine Upwind");
  sprintf(key, "Solver.TerrainFollowingGrid.SlopeUpwindFormulation");
  switch_name = GetStringDefault(key, "Original");
  switch_value = NA_NameToIndex(upwind_switch_na, switch_name);
  switch (switch_value)
  {
    case 0:
    {
      public_xtra->tfgupwind = 0;
      break;
    }

    case 1:
    {
      public_xtra->tfgupwind = 1;
      break;
    }

    case 2:
    {
      public_xtra->tfgupwind = 2;
      break;
    }

    default:
    {
      InputError("Error: Invalid value <%s> for key <%s>\n", switch_name,
                 key);
    }
  }
  NA_FreeNameArray(upwind_switch_na);

  switch_na = NA_NewNameArray("False True");
  sprintf(key, "Solver.Nonlinear.UseJacobian");
  switch_name = GetStringDefault(key, "False");
  switch_value = NA_NameToIndex(switch_na, switch_name);
  switch (switch_value)
  {
    case 0:
    {
      public_xtra->type = no_nonlinear_jacobian;
      break;
    }

    case 1:
    {
      public_xtra->type = not_set;
      break;
    }

    default:
    {
      InputError("Error: Invalid value <%s> for key <%s>\n", switch_name,
                 key);
    }
  }
  NA_FreeNameArray(switch_na);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * RichardsJacobianEvalFreePublicXtra
 *--------------------------------------------------------------------------*/

void  RichardsJacobianEvalFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);


  if (public_xtra)
  {
    tfree(public_xtra);
  }
}


/*--------------------------------------------------------------------------
 * RichardsJacobianEvalSizeOfTempData
 *--------------------------------------------------------------------------*/

int  RichardsJacobianEvalSizeOfTempData()
{
  int sz = 0;

  return sz;
}
