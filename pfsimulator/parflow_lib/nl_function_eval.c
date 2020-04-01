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
#include "parflow_config.h"

#ifdef USING_PARALLEL
extern "C"{
#endif

#include "parflow.h"
#include "pf_parallel.h"
#include "llnlmath.h"
#include "llnltyps.h"
//#include "math.h"
#include "float.h"


/*---------------------------------------------------------------------
 * Define module structures
 *---------------------------------------------------------------------*/

typedef struct {
  int time_index;
  double SpinupDampP1;      // NBE
  double SpinupDampP2;      // NBE
  int tfgupwind;           //@RMM added for TFG formulation switch
} PublicXtra;

typedef struct {
  Problem      *problem;

  PFModule     *density_module;
  PFModule     *saturation_module;
  PFModule     *rel_perm_module;
  PFModule     *phase_source;
  PFModule     *bc_pressure;
  PFModule     *bc_internal;
  PFModule     *overlandflow_module;  //DOK
  PFModule     *overlandflow_module_diff;  //@RMM
  PFModule     *overlandflow_module_kin;
} InstanceXtra;

/*---------------------------------------------------------------------
 * Define macros for function evaluation
 *---------------------------------------------------------------------*/
#ifdef HAVE_OMP

#define PMean(a, b, c, d)   _HarmonicMean(c, d)
#define PMeanDZ(a, b, c, d) _HarmonicMeanDZ(a, b, c, d)
#define RPMean(a, b, c, d)  _UpstreamMean(a, b, c, d)
#define Mean(a, b)          _ArithmeticMean(a, b)

#else

#define PMean(a, b, c, d)    HarmonicMean(c, d)
#define PMeanDZ(a, b, c, d)     HarmonicMeanDZ(a, b, c, d)
#define RPMean(a, b, c, d)   UpstreamMean(a, b, c, d)
#define Mean(a, b)            ArithmeticMean(a, b)

#endif

/*  This routine provides the interface between KINSOL and ParFlow
 *  for function evaluations.  */

void     KINSolFunctionEval(
                            int      size,
                            N_Vector pressure,
                            N_Vector fval,
                            void *   current_state)
{
  PFModule  *nl_function_eval = StateFunc(((State*)current_state));
  ProblemData *problem_data = StateProblemData(((State*)current_state));
  Vector      *old_pressure = StateOldPressure(((State*)current_state));
  Vector      *saturation = StateSaturation(((State*)current_state));
  Vector      *old_saturation = StateOldSaturation(((State*)current_state));
  Vector      *density = StateDensity(((State*)current_state));
  Vector      *old_density = StateOldDensity(((State*)current_state));
  double dt = StateDt(((State*)current_state));
  double time = StateTime(((State*)current_state));
  Vector       *evap_trans = StateEvapTrans(((State*)current_state));
  Vector       *ovrl_bc_flx = StateOvrlBcFlx(((State*)current_state));

  /* velocity vectors jjb */
  Vector       *x_velocity = StateXvel(((State*)current_state));
  Vector       *y_velocity = StateYvel(((State*)current_state));
  Vector       *z_velocity = StateZvel(((State*)current_state));

  (void)size;

  PFModuleInvokeType(NlFunctionEvalInvoke, nl_function_eval,
                     (pressure, fval, problem_data, saturation, old_saturation,
                      density, old_density, dt, time, old_pressure, evap_trans,
                      ovrl_bc_flx, x_velocity, y_velocity, z_velocity));

  return;
}


/*  This routine evaluates the nonlinear function based on the current
 *  pressure values.  This evaluation is basically an application
 *  of the stencil to the pressure array. */

void NlFunctionEval(Vector *     pressure, /* Current pressure values */
                    Vector *     fval, /* Return values of the nonlinear function */
                    ProblemData *problem_data,  /* Geometry data for problem */
                    Vector *     saturation, /* Saturation / work vector */
                    Vector *     old_saturation, /* Saturation values at previous time step */
                    Vector *     density, /* Density vector */
                    Vector *     old_density, /* Density values at previous time step */
                    double       dt, /* Time step size */
                    double       time, /* New time value */
                    Vector *     old_pressure,
                    Vector *     evap_trans, /*sk sink term from land surface model*/
                    Vector *     ovrl_bc_flx, /*sk overland flow boundary fluxes*/
                    Vector *     x_velocity, /* velocity vectors jjb */
                    Vector *     y_velocity,
                    Vector *     z_velocity)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  Problem     *problem = (instance_xtra->problem);

  PFModule    *density_module = (instance_xtra->density_module);
  PFModule    *saturation_module = (instance_xtra->saturation_module);
  PFModule    *rel_perm_module = (instance_xtra->rel_perm_module);
  PFModule    *phase_source = (instance_xtra->phase_source);
  PFModule    *bc_pressure = (instance_xtra->bc_pressure);
  PFModule    *bc_internal = (instance_xtra->bc_internal);
  PFModule    *overlandflow_module = (instance_xtra->overlandflow_module);
  PFModule    *overlandflow_module_diff = (instance_xtra->overlandflow_module_diff);
  PFModule    *overlandflow_module_kin = (instance_xtra->overlandflow_module_kin);

  Vector      *porosity = ProblemDataPorosity(problem_data);
  Vector      *permeability_x = ProblemDataPermeabilityX(problem_data);
  Vector      *permeability_y = ProblemDataPermeabilityY(problem_data);
  Vector      *permeability_z = ProblemDataPermeabilityZ(problem_data);
  Vector      *sstorage = ProblemDataSpecificStorage(problem_data);            //sk
  Vector      *x_sl = ProblemDataTSlopeX(problem_data);                //sk
  Vector      *y_sl = ProblemDataTSlopeY(problem_data);                //sk
  Vector      *man = ProblemDataMannings(problem_data);                 //sk

  /* @RMM terrain following grid slope variables */
  Vector      *x_ssl = ProblemDataSSlopeX(problem_data);               //@RMM
  Vector      *y_ssl = ProblemDataSSlopeY(problem_data);               //@RMM

  /* @RMM variable dz multiplier */
  Vector      *z_mult = ProblemDataZmult(problem_data);              //@RMM

/* @RMM Flow Barrier / Boundary values */
  Vector      *FBx = ProblemDataFBx(problem_data);
  Vector      *FBy = ProblemDataFBy(problem_data);
  Vector      *FBz = ProblemDataFBz(problem_data);

  double gravity = ProblemGravity(problem);
  double viscosity = ProblemPhaseViscosity(problem, 0);

  Grid        *grid = VectorGrid(pressure);
  Grid        *grid2d = VectorGrid(x_sl);
  GrGeomSolid *gr_domain = ProblemDataGrDomain(problem_data);

  BeginTiming(public_xtra->time_index);

  /* diffusive test here, this is NOT PF style and should be
   * re-done putting keys in BC Pressure Package and adding to the
   * datastructure for overlandflowBC */
  int diffusive;             //@RMM
  diffusive = GetIntDefault("OverlandFlowDiffusive", 0);

  int overlandspinup;              //@RMM
  overlandspinup = GetIntDefault("OverlandFlowSpinUp", 0);

  /* Pass pressure values to neighbors.  */
  VectorUpdateCommHandle  *handle;
  handle = InitVectorUpdate(pressure, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  Vector      *KW, *KE, *KN, *KS;
  Vector      *qx, *qy;
  KW = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);
  KE = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);
  KN = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);
  KS = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);
  qx = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);
  qy = NewVectorType(grid2d, 1, 1, vector_cell_centered_2D);


/* Initialize function values to zero. */
  PFVConstInit(0.0, fval);

  BCStruct    *bc_struct;
  bc_struct = PFModuleInvokeType(BCPressureInvoke, bc_pressure,
                                 (problem_data, grid, gr_domain, time));

  Vector *u_right_vec, *u_front_vec, *u_upper_vec; // @IJB
  u_right_vec = NewVectorType(grid, 1, 1, fval->type);
  u_front_vec = NewVectorType(grid, 1, 1, fval->type);
  u_upper_vec = NewVectorType(grid, 1, 1, fval->type);

  /* Re-use saturation vector to save memory */
  Vector      *rel_perm = saturation;
  Vector      *source = saturation;

#pragma omp parallel
  {

    Subvector   *u_right_sub, *u_front_sub, *u_upper_sub; // @IJB
    double      *u_right_dat, *u_front_dat, *u_upper_dat; // @IJB

  /* Overland flow variables */  //sk
  Subvector   *kw_sub, *ke_sub, *kn_sub, *ks_sub, *qx_sub, *qy_sub;
  Subvector   *x_sl_sub, *y_sl_sub, *mann_sub;
  Subvector   *obf_sub;
  double      *kw_, *ke_, *kn_, *ks_, *qx_, *qy_;
  double      *x_sl_dat, *y_sl_dat, *mann_dat;
  double      *obf_dat;
  double q_overlnd;
  double sep;          // scaling difference temp var @RMM

  Subvector   *x_ssl_sub, *y_ssl_sub;    //@RMM
  double      *x_ssl_dat, *y_ssl_dat;     //@RMM

  Subvector   *z_mult_sub;    //@RMM
  double      *z_mult_dat;    //@RMM


  Subvector   *FBx_sub, *FBy_sub, *FBz_sub;  //@RMM
  double      *FBx_dat, *FBy_dat, *FBz_dat;   //@RMM

  Subgrid     *subgrid;

  Subvector   *p_sub, *d_sub, *od_sub, *s_sub, *os_sub, *po_sub, *op_sub, *ss_sub, *et_sub;
  Subvector   *f_sub, *rp_sub, *permx_sub, *permy_sub, *permz_sub;

  Subvector   *vx_sub, *vy_sub, *vz_sub;  //jjb
  double      *vx, *vy, *vz;  //jjb
  int vxi, vyi, vzi;         //jjb



  double      *pp, *odp, *sp, *osp, *pop, *fp, *dp, *rpp, *opp, *ss, *et;
  double      *permxp, *permyp, *permzp;

  int i, j, k, r, is;
  int ix, iy, iz;
  int nx, ny, nz, gnx, gny;
  int nx_f, ny_f, nz_f;
  int nx_p, ny_p, nz_p;
  int nx_po, ny_po, nz_po;
  int sy_p, sz_p;
  int ip, ipo, io;


  double dx, dy, dz, vol, ffx, ffy, ffz;
  double u_right, u_front, u_upper;
  double diff = 0.0e0;
  double updir = 0.0e0;
  double lower_cond, upper_cond;
  //@RMM : terms for gravity/terrain
  double x_dir_g, y_dir_g, z_dir_g, del_x_slope, del_y_slope, x_dir_g_c, y_dir_g_c;

  double      *bc_patch_values;
  double u_old = 0.0e0;
  double u_new = 0.0e0;
  double value;
  int         *fdir;
  int ipatch, ival;
  int dir = 0;


  /* Calculate pressure dependent properties: density and saturation */
  double dtmp;
  PFModuleInvokeType(PhaseDensityInvoke, density_module, (0, pressure, density, &dtmp, &dtmp,
                                                          CALCFCN));
  /* @MCB:
     This barrier is necessary because the PhaseDensity module has no barriers.
     RichardsJacobianEval calls this module twice, for CALCFCN and CALCDER, and so
     the lack of barriers in the module is beneficial there.
     This is no more expensive than having the barrier inside the module in this case.
  */
  BARRIER;

  PFModuleInvokeType(SaturationInvoke, saturation_module, (saturation, pressure, density,
                                                           gravity, problem_data, CALCFCN));

  /* Calculate accumulation terms for the function values */

  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);
    Subgrid* grid2d_subgrid = GridSubgrid(grid2d, is);
    int grid2d_iz = SubgridIZ(grid2d_subgrid);

    d_sub = VectorSubvector(density, is);
    od_sub = VectorSubvector(old_density, is);
    p_sub = VectorSubvector(pressure, is);
    op_sub = VectorSubvector(old_pressure, is);
    s_sub = VectorSubvector(saturation, is);
    os_sub = VectorSubvector(old_saturation, is);
    po_sub = VectorSubvector(porosity, is);
    f_sub = VectorSubvector(fval, is);

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

    /* @RMM added to provide access FB values */
    FBx_sub = VectorSubvector(FBx, is);
    FBy_sub = VectorSubvector(FBy, is);
    FBz_sub = VectorSubvector(FBz, is);

    /* @RMM added to provide FB values */
    FBx_dat = SubvectorData(FBx_sub);
    FBy_dat = SubvectorData(FBy_sub);
    FBz_dat = SubvectorData(FBz_sub);

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

    nx_f = SubvectorNX(f_sub);
    ny_f = SubvectorNY(f_sub);
    nz_f = SubvectorNZ(f_sub);

    nx_po = SubvectorNX(po_sub);
    ny_po = SubvectorNY(po_sub);
    nz_po = SubvectorNZ(po_sub);

    dp = SubvectorData(d_sub);
    odp = SubvectorData(od_sub);
    sp = SubvectorData(s_sub);
    pp = SubvectorData(p_sub);
    opp = SubvectorData(op_sub);
    osp = SubvectorData(os_sub);
    pop = SubvectorData(po_sub);
    fp = SubvectorData(f_sub);

    _GrGeomInLoop(NoWait, NO_LOCALS,
                  i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      ip = SubvectorEltIndex(f_sub, i, j, k);
      ipo = SubvectorEltIndex(po_sub, i, j, k);
      io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

      /*     del_x_slope = (1.0/cos(atan(x_ssl_dat[io])));
       *   del_y_slope = (1.0/cos(atan(y_ssl_dat[io])));  */
      del_x_slope = 1.0;
      del_y_slope = 1.0;

      fp[ip] = (sp[ip] * dp[ip] - osp[ip] * odp[ip]) * pop[ipo] * vol * del_x_slope * del_y_slope * z_mult_dat[ip];
    });
  }

  /*@ Add in contributions from compressible storage */
  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);
    Subgrid       *grid2d_subgrid = GridSubgrid(grid2d, is);
    int grid2d_iz = SubgridIZ(grid2d_subgrid);

    ss_sub = VectorSubvector(sstorage, is);

    d_sub = VectorSubvector(density, is);
    od_sub = VectorSubvector(old_density, is);
    p_sub = VectorSubvector(pressure, is);
    op_sub = VectorSubvector(old_pressure, is);
    s_sub = VectorSubvector(saturation, is);
    os_sub = VectorSubvector(old_saturation, is);
    f_sub = VectorSubvector(fval, is);

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

    nx_f = SubvectorNX(f_sub);
    ny_f = SubvectorNY(f_sub);
    nz_f = SubvectorNZ(f_sub);

    ss = SubvectorData(ss_sub);

    dp = SubvectorData(d_sub);
    odp = SubvectorData(od_sub);
    sp = SubvectorData(s_sub);
    pp = SubvectorData(p_sub);
    opp = SubvectorData(op_sub);
    osp = SubvectorData(os_sub);
    fp = SubvectorData(f_sub);

    _GrGeomInLoop(InParallel, NO_LOCALS,
                  i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      ip = SubvectorEltIndex(f_sub, i, j, k);
      io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

      /*   del_x_slope = (1.0/cos(atan(x_ssl_dat[io])));
       * del_y_slope = (1.0/cos(atan(y_ssl_dat[io])));  */
      del_x_slope = 1.0;
      del_y_slope = 1.0;
      fp[ip] += ss[ip] * vol * del_x_slope * del_y_slope * z_mult_dat[ip] *
                (pp[ip] * sp[ip] * dp[ip] - opp[ip] * osp[ip] * odp[ip]);
    });
  }

  /* Add in contributions from source terms - user specified sources and
   * flux wells.  Calculate phase source values overwriting current
   * saturation vector */
  /* @MCB:
     This module has implicit barriers since it is only called once
  */
  PFModuleInvokeType(PhaseSourceInvoke, phase_source, (source, 0, problem, problem_data,
                                                       time));

  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);
    Subgrid       *grid2d_subgrid = GridSubgrid(grid2d, is);
    int grid2d_iz = SubgridIZ(grid2d_subgrid);

    s_sub = VectorSubvector(source, is);
    f_sub = VectorSubvector(fval, is);
    et_sub = VectorSubvector(evap_trans, is);

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

    nx_f = SubvectorNX(f_sub);
    ny_f = SubvectorNY(f_sub);
    nz_f = SubvectorNZ(f_sub);

    sp = SubvectorData(s_sub);
    fp = SubvectorData(f_sub);
    et = SubvectorData(et_sub);

    /* @RMM added to provide access to x/y slopes */
    x_ssl_sub = VectorSubvector(x_ssl, is);
    y_ssl_sub = VectorSubvector(y_ssl, is);
    /* @RMM  added to provide slopes to terrain fns */
    x_ssl_dat = SubvectorData(x_ssl_sub);
    y_ssl_dat = SubvectorData(y_ssl_sub);
    /* @RMM added to provide access to zmult */
    z_mult_sub = VectorSubvector(z_mult, is);
    /* @RMM added to provide variable dz */
    z_mult_dat = SubvectorData(z_mult_sub);
    /* @RMM added to provide access FB values */
    FBx_sub = VectorSubvector(FBx, is);
    FBy_sub = VectorSubvector(FBy, is);
    FBz_sub = VectorSubvector(FBz, is);

    /* @RMM added to provide FB values */
    FBx_dat = SubvectorData(FBx_sub);
    FBy_dat = SubvectorData(FBy_sub);
    FBz_dat = SubvectorData(FBz_sub);

    /* @MCB:
       This doesn't need to wait because the next loop (BCStructPatch) is only writing to pressure
       but we aren't reading from pressure here
    */
    _GrGeomInLoop(NoWait, NO_LOCALS,
                  i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      ip = SubvectorEltIndex(f_sub, i, j, k);
      io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

      /* del_x_slope = (1.0/cos(atan(x_ssl_dat[io])));
       * del_y_slope = (1.0/cos(atan(y_ssl_dat[io])));  */
      del_x_slope = 1.0;
      del_y_slope = 1.0;
      fp[ip] -= vol * del_x_slope * del_y_slope * z_mult_dat[ip] * dt * (sp[ip] + et[ip]);
    });
  }


  /*
   * Temporarily insert boundary pressure values for Dirichlet
   * boundaries into cells that are in the inactive region but next
   * to a Dirichlet boundary condition.  These values are required
   * for use in the rel_perm_module to compute rel_perm values for
   * these cells. They needed for upstream weighting in mobilities.
   *
   * NOTES:
   *
   * These values must be later removed from the pressure field and
   * fval needs to be adjusted for these cells to make the inactive
   * region problem decoupled from the active region cells for the
   * solver.
   *
   * Densities are currently defined everywhere so should be valid for
   * these boundary cells.
   *
   * SGS not sure if this will work or not so left it here for later
   * exploration.  This is a little hacky in the sense that we are
   * inserting values and then need to overwrite them again.  It
   * might be more clean to rewrite the Dirichlet boundary condition
   * code to not require the values be in the pressure field for
   * these cells but instead grab the values out of the
   * BCStructPatchValues as was done here.  In other words use
   * bc_patch_values[ival] in rel_perm_module code and remove this
   * loop.
   */


  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    p_sub = VectorSubvector(pressure, is);

    nx_p = SubvectorNX(p_sub);
    ny_p = SubvectorNY(p_sub);
    nz_p = SubvectorNZ(p_sub);

    sy_p = nx_p;
    sz_p = ny_p * nx_p;

    pp = SubvectorData(p_sub);

    int pp_idx = 0;
    ForBCStructNumPatches(ipatch, bc_struct)
    {
      bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);
      ForPatchCellsPerFace(DirichletBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup({
                               pp_idx = 0;
                               ip = SubvectorEltIndex(p_sub, i, j, k);
                               value = bc_patch_values[ival];
                             }),
                           FACE(Left,    { pp_idx = ip - 1; }),
                           FACE(Right,   { pp_idx = ip + 1; }),
                           FACE(Down,    { pp_idx = ip - sy_p; }),
                           FACE(Up,      { pp_idx = ip + sy_p; }),
                           FACE(Back,    { pp_idx = ip - sz_p; }),
                           FACE(Front,   { pp_idx = ip + sz_p; }),
                           CellFinalize( { pp[pp_idx] = value; }),
                           AfterAllCells(DoNothing)
        );
    }          /* End ipatch loop */
  }            /* End subgrid loop */

  /* Calculate relative permeability values overwriting current
   * phase source values */

  PFModuleInvokeType(PhaseRelPermInvoke, rel_perm_module,
                     (rel_perm, pressure, density, gravity, problem_data,
                      CALCFCN));
  BARRIER;

  /* Calculate contributions from second order derivatives and gravity */
  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    /* velocity vectors jjb */
    vx_sub = VectorSubvector(x_velocity, is);
    vy_sub = VectorSubvector(y_velocity, is);
    vz_sub = VectorSubvector(z_velocity, is);

    Subgrid       *grid2d_subgrid = GridSubgrid(grid2d, is);
    int grid2d_iz = SubgridIZ(grid2d_subgrid);

    p_sub = VectorSubvector(pressure, is);
    d_sub = VectorSubvector(density, is);
    rp_sub = VectorSubvector(rel_perm, is);
    f_sub = VectorSubvector(fval, is);
    permx_sub = VectorSubvector(permeability_x, is);
    permy_sub = VectorSubvector(permeability_y, is);
    permz_sub = VectorSubvector(permeability_z, is);
    /* @RMM added to provide access to x/y slopes */
    x_ssl_sub = VectorSubvector(x_ssl, is);
    y_ssl_sub = VectorSubvector(y_ssl, is);

    /* @RMM added to provide access to zmult */
    z_mult_sub = VectorSubvector(z_mult, is);


    /* @IJB added to enable parallelism */
    u_right_sub = VectorSubvector(u_right_vec, is);
    u_front_sub = VectorSubvector(u_front_vec, is);
    u_upper_sub = VectorSubvector(u_upper_vec, is);

    /* RDF: assumes resolutions are the same in all 3 directions */
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

    nx_p = SubvectorNX(p_sub);
    ny_p = SubvectorNY(p_sub);
    nz_p = SubvectorNZ(p_sub);

    int sx_p = 1;
    sy_p = nx_p;
    sz_p = ny_p * nx_p;

    /* velocity accessors jjb */
    vx = SubvectorData(vx_sub);
    vy = SubvectorData(vy_sub);
    vz = SubvectorData(vz_sub);

    pp = SubvectorData(p_sub);
    dp = SubvectorData(d_sub);
    rpp = SubvectorData(rp_sub);
    fp = SubvectorData(f_sub);
    permxp = SubvectorData(permx_sub);
    permyp = SubvectorData(permy_sub);
    permzp = SubvectorData(permz_sub);

    /* @RMM  added to provide slopes to terrain fns */
    x_ssl_dat = SubvectorData(x_ssl_sub);
    y_ssl_dat = SubvectorData(y_ssl_sub);

    /* @RMM added to provide variable dz */
    z_mult_dat = SubvectorData(z_mult_sub);


    /* @IJB added to enable parallelism */
    u_right_dat = SubvectorData(u_right_sub);
    u_front_dat = SubvectorData(u_front_sub);
    u_upper_dat = SubvectorData(u_upper_sub);

    qx_sub = VectorSubvector(qx, is);

    _GrGeomInLoop(InParallel, NO_LOCALS,
                  i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      ip = SubvectorEltIndex(p_sub, i, j, k);
      io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

      /* @RMM: modified the terrain-following transform
       * to be swtichable in the UZ
       * terms:
       * 1. x dir terrain tendency:  gravity*sin(atan(x_ssl_dat[io]))
       * 2. y dir terrain tendency:  gravity*sin(atan(y_ssl_dat[io]))
       * 3. change in delta-x due to slope: (1.0/cos(atan(x_ssl_dat[io])))
       * 4. change in delta-y due to slope: (1.0/cos(atan(y_ssl_dat[io])))
       * Depending on formulation chosen slopes are either assumed to be cell-centered or
       * upwind
       */

      /* velocity subvector indices jjb */
      vxi = SubvectorEltIndex(vx_sub, i + 1, j, k);
      vyi = SubvectorEltIndex(vy_sub, i, j + 1, k);
      vzi = SubvectorEltIndex(vz_sub, i, j, k + 1);

      z_dir_g = 1.0;
      del_x_slope = 1.0;
      del_y_slope = 1.0;

//@RMM  tfgupwind == 0 (default) should give original behavior
// tfgupwind 1 should still use sine but upwind
// tfgupwdin 2 just upwind
      switch (public_xtra->tfgupwind)
      {
        case 0:
          {
            // default formulation in Maxwell 2013
            x_dir_g = Mean(gravity * sin(atan(x_ssl_dat[io])), gravity * sin(atan(x_ssl_dat[io + 1])));
            x_dir_g_c = Mean(gravity * cos(atan(x_ssl_dat[io])), gravity * cos(atan(x_ssl_dat[io + 1])));
            y_dir_g = Mean(gravity * sin(atan(y_ssl_dat[io])), gravity * sin(atan(y_ssl_dat[io + sy_p])));
            y_dir_g_c = Mean(gravity * cos(atan(y_ssl_dat[io])), gravity * cos(atan(y_ssl_dat[io + sy_p])));
            break;
          }

        case 1:
          {
            // direct upwinding, no averaging with sines
            x_dir_g = gravity * sin(atan(x_ssl_dat[io]));
            x_dir_g_c = gravity * cos(atan(x_ssl_dat[io]));
            y_dir_g = gravity * sin(atan(y_ssl_dat[io]));
            y_dir_g_c = gravity * cos(atan(y_ssl_dat[io]));
            break;
          }

        case 2:
          {
            // direct upwinding, no averaging no sines
            x_dir_g = x_ssl_dat[io];
            x_dir_g_c = 1.0;
            y_dir_g = y_ssl_dat[io];
            y_dir_g_c = 1.0;
            break;
          }
      }
      /* Calculate right face velocity.
       * diff >= 0 implies flow goes left to right */

      diff = pp[ip] - pp[ip + 1];
      updir = (diff / dx) * x_dir_g_c - x_dir_g;

      u_right = z_mult_dat[ip] * ffx * del_y_slope * PMean(pp[ip], pp[ip + 1],
                                                           permxp[ip], permxp[ip + 1])
                * (diff / (dx * del_x_slope)) * x_dir_g_c
                * RPMean(updir, 0.0,
                         rpp[ip] * dp[ip],
                         rpp[ip + 1] * dp[ip + 1])
                / viscosity;

      /* Calculate right face velocity gravity terms
       * @RMM added sin* g term to test terrain-following grid
       * upwind on pressure is currently implemented
       * Sx < 0 implies flow goes left to right */

      u_right += z_mult_dat[ip] * ffx * del_y_slope * PMean(pp[ip], pp[ip + 1],
                                                            permxp[ip], permxp[ip + 1])
                 * (-x_dir_g)
                 * RPMean(updir, 0.0, rpp[ip] * dp[ip],
                          rpp[ip + 1] * dp[ip + 1])
                 / viscosity;


      /* Calculate front face velocity.
       * diff >= 0 implies flow goes back to front */
      diff = pp[ip] - pp[ip + sy_p];
      updir = (diff / dy) * y_dir_g_c - y_dir_g;

      u_front = z_mult_dat[ip] * ffy * del_x_slope
                * PMean(pp[ip], pp[ip + sy_p], permyp[ip], permyp[ip + sy_p])
                * (diff / (dy * del_y_slope)) * y_dir_g_c
                * RPMean(updir, 0.0,
                         rpp[ip] * dp[ip],
                         rpp[ip + sy_p] * dp[ip + sy_p])
                / viscosity;

      /* Calculate front face velocity gravity terms
       * @RMM added sin* g term to test terrain-following grid
       * note upwinding on gravity terms not pressure
       * Sy < 0 implies flow goes from left to right
       */

      u_front += z_mult_dat[ip] * ffy * del_x_slope
                 * PMean(pp[ip], pp[ip + sy_p], permyp[ip], permyp[ip + sy_p])
                 * (-y_dir_g)
                 * RPMean(updir, 0.0, rpp[ip] * dp[ip],
                          rpp[ip + sy_p] * dp[ip + sy_p])
                 / viscosity;

      /* Calculate upper face velocity.
       * diff >= 0 implies flow goes lower to upper
       */
      sep = dz * (Mean(z_mult_dat[ip], z_mult_dat[ip + sz_p]));


      lower_cond = pp[ip] / sep
                   - (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p]))
                   * dp[ip] * gravity * z_dir_g;

      upper_cond = pp[ip + sz_p] / sep
                   + (z_mult_dat[ip + sz_p] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p]))
                   * dp[ip + sz_p] * gravity * z_dir_g;


      diff = (lower_cond - upper_cond);

      u_upper = ffz * del_x_slope * del_y_slope
                * PMeanDZ(permzp[ip], permzp[ip + sz_p], z_mult_dat[ip], z_mult_dat[ip + sz_p])
                * diff
                * RPMean(lower_cond, upper_cond, rpp[ip] * dp[ip],
                         rpp[ip + sz_p] * dp[ip + sz_p])
                / viscosity;

/*  add in flow barrier values
 * assumes that ip is the cell face between ip and ip+1
 * ip and ip+sy_p and ip and ip+sz_p in x, y, z directions */
      u_right = u_right * FBx_dat[ip];
      u_front = u_front * FBy_dat[ip];
      u_upper = u_upper * FBz_dat[ip];

      /* velocity data jjb */
      vx[vxi] = u_right / ffx;
      vy[vyi] = u_front / ffy;
      vz[vzi] = u_upper / ffz;

      u_right_dat[ip] = u_right;
      u_front_dat[ip] = u_front;
      u_upper_dat[ip] = u_upper;

      /*
      fp[ip] += dt * (u_right + u_front + u_upper);
      fp[ip + 1] -= dt * u_right;
      fp[ip + sy_p] -= dt * u_front;
      fp[ip + sz_p] -= dt * u_upper;
      */
      /*
      PlusEquals(fp[ip], (dt * (u_right + u_front + u_upper)));
      PlusEquals(fp[ip + 1], -(dt * u_right));
      PlusEquals(fp[ip + sy_p], -(dt * u_front));
      PlusEquals(fp[ip + sz_p], -(dt * u_upper));
      */
    });
  }

  #pragma omp single
  {
// Gather portion
  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    p_sub = VectorSubvector(pressure, is);
    f_sub = VectorSubvector(fval, is);

    /* @IJB added to enable parallelism */
    u_right_sub = VectorSubvector(u_right_vec, is);
    u_front_sub = VectorSubvector(u_front_vec, is);
    u_upper_sub = VectorSubvector(u_upper_vec, is);

    ix = SubgridIX(subgrid) - 1;
    iy = SubgridIY(subgrid) - 1;
    iz = SubgridIZ(subgrid) - 1;

    nx = SubgridNX(subgrid) + 1;
    ny = SubgridNY(subgrid) + 1;
    nz = SubgridNZ(subgrid) + 1;

    nx_p = SubvectorNX(p_sub);
    ny_p = SubvectorNY(p_sub);
    nz_p = SubvectorNZ(p_sub);

    int sx_p = 1;
    sy_p = nx_p;
    sz_p = ny_p * nx_p;

    fp = SubvectorData(f_sub);

    /* @IJB added to enable parallelism */
    u_right_dat = SubvectorData(u_right_sub);
    u_front_dat = SubvectorData(u_front_sub);
    u_upper_dat = SubvectorData(u_upper_sub);

    //_GrGeomInLoop(InParallel, NO_LOCALS,
    GrGeomInLoop(
                  i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,
    {
      ip = SubvectorEltIndex(p_sub, i, j, k);
      fp[ip] += dt * (u_right_dat[ip] + u_front_dat[ip] + u_upper_dat[ip]);
      fp[ip + 1] -= dt * u_right_dat[ip];
      fp[ip + sy_p] -= dt * u_front_dat[ip];
      fp[ip + sz_p] -= dt * u_upper_dat[ip];
    });
  }
  }

  /*  Calculate correction for boundary conditions */


  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);
    Subgrid       *grid2d_subgrid = GridSubgrid(grid2d, is);
    int grid2d_iz = SubgridIZ(grid2d_subgrid);

    d_sub = VectorSubvector(density, is);
    rp_sub = VectorSubvector(rel_perm, is);
    f_sub = VectorSubvector(fval, is);
    permx_sub = VectorSubvector(permeability_x, is);
    permy_sub = VectorSubvector(permeability_y, is);
    permz_sub = VectorSubvector(permeability_z, is);

    p_sub = VectorSubvector(pressure, is);
    op_sub = VectorSubvector(old_pressure, is);
    os_sub = VectorSubvector(old_saturation, is);

    /* @RMM added to provide access to x/y slopes */
    x_ssl_sub = VectorSubvector(x_ssl, is);
    y_ssl_sub = VectorSubvector(y_ssl, is);

    // sk Overland flow
    kw_sub = VectorSubvector(KW, is);
    ke_sub = VectorSubvector(KE, is);
    kn_sub = VectorSubvector(KN, is);
    ks_sub = VectorSubvector(KS, is);
    qx_sub = VectorSubvector(qx, is);
    qy_sub = VectorSubvector(qy, is);
    x_sl_sub = VectorSubvector(x_sl, is);
    y_sl_sub = VectorSubvector(y_sl, is);
    mann_sub = VectorSubvector(man, is);
    /*
     * SGS TODO This looks very wrong, why going to DB here, should
     * come from DS
     */
    gnx = GetInt("ComputationalGrid.NX");
    gny = GetInt("ComputationalGrid.NY");
    obf_sub = VectorSubvector(ovrl_bc_flx, is);


    dx = SubgridDX(subgrid);
    dy = SubgridDY(subgrid);
    dz = SubgridDZ(subgrid);

    nx = SubgridNX(subgrid);
    ny = SubgridNY(subgrid);

    ix = SubgridIX(subgrid);
    iy = SubgridIY(subgrid);

    ffx = dy * dz;
    ffy = dx * dz;
    ffz = dx * dy;

    vol = dx * dy * dz;

    nx_p = SubvectorNX(p_sub);
    ny_p = SubvectorNY(p_sub);
    nz_p = SubvectorNZ(p_sub);

    sy_p = nx_p;
    sz_p = ny_p * nx_p;

    dp = SubvectorData(d_sub);
    rpp = SubvectorData(rp_sub);
    fp = SubvectorData(f_sub);
    permxp = SubvectorData(permx_sub);
    permyp = SubvectorData(permy_sub);
    permzp = SubvectorData(permz_sub);

    kw_ = SubvectorData(kw_sub);
    ke_ = SubvectorData(ke_sub);
    kn_ = SubvectorData(kn_sub);
    ks_ = SubvectorData(ks_sub);
    qx_ = SubvectorData(qx_sub);
    qy_ = SubvectorData(qy_sub);
    x_sl_dat = SubvectorData(x_sl_sub);
    y_sl_dat = SubvectorData(y_sl_sub);
    mann_dat = SubvectorData(mann_sub);
    obf_dat = SubvectorData(obf_sub);

    pp = SubvectorData(p_sub);
    opp = SubvectorData(op_sub);
    osp = SubvectorData(os_sub);

    /* @RMM  added to provide slopes to terrain fns */
    x_ssl_dat = SubvectorData(x_ssl_sub);
    y_ssl_dat = SubvectorData(y_ssl_sub);
    /* @RMM added to provide access to zmult */
    z_mult_sub = VectorSubvector(z_mult, is);
    /* @RMM added to provide variable dz */
    z_mult_dat = SubvectorData(z_mult_sub);



    ForBCStructNumPatches(ipatch, bc_struct)
    {
      bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);

      ForPatchCellsPerFace(DirichletBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup(
                           {
                             ip = SubvectorEltIndex(p_sub, i, j, k);
                             io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

                             value = bc_patch_values[ival];
                             x_dir_g = 0.0;
                             y_dir_g = 0.0;
                             z_dir_g = 1.0;

                             del_x_slope = (1.0 / cos(atan(x_ssl_dat[io])));
                             del_y_slope = (1.0 / cos(atan(y_ssl_dat[io])));

                             del_x_slope = 1.0;
                             del_y_slope = 1.0;
                           }),
                           FACE(Left,
                           {
                             dir = -1;
                             diff = pp[ip - 1] - pp[ip];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip - 1], pp[ip], permxp[ip - 1], permxp[ip])
                                     * (diff / dx * del_x_slope)
                                     * RPMean(pp[ip - 1], pp[ip],
                                              rpp[ip - 1] * dp[ip - 1], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope *
                                      PMean(pp[ip - 1], pp[ip],
                                            permxp[ip - 1], permxp[ip])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip - 1], pp[ip], rpp[ip - 1] * dp[ip - 1],
                                               rpp[ip] * dp[ip])
                                      / viscosity;

                             diff = value - pp[ip];
                             u_new = RPMean(value, pp[ip],
                                            rpp[ip - 1] * dp[ip - 1], rpp[ip] * dp[ip]);
                             u_new = u_new * z_mult_dat[ip] * ffx * del_y_slope
                                     * (permxp[ip] / viscosity)
                                     * 2.0 * (diff / dx);
                           }),
                           FACE(Right,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + 1];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip], pp[ip + 1], permxp[ip], permxp[ip + 1])
                                     * (diff / dx * del_x_slope)
                                     * RPMean(pp[ip], pp[ip + 1],
                                              rpp[ip] * dp[ip], rpp[ip + 1] * dp[ip + 1])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip], pp[ip + 1],
                                              permxp[ip], permxp[ip + 1])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip], pp[ip + 1], rpp[ip] * dp[ip],
                                               rpp[ip + 1] * dp[ip + 1])
                                      / viscosity;

                             diff = pp[ip] - value;
                             u_new = RPMean(pp[ip], value,
                                            rpp[ip] * dp[ip], rpp[ip + 1] * dp[ip + 1]);
                             u_new = u_new * z_mult_dat[ip] * ffx * del_y_slope
                                     * (permxp[ip] / viscosity)
                                     * 2.0 * (diff / dx);
                           }),
                           FACE(Down,
                           {
                             dir = -1;
                             diff = pp[ip - sy_p] - pp[ip];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip - sy_p], pp[ip],
                                             permyp[ip - sy_p], permyp[ip])
                                     * (diff / dy * del_y_slope)
                                     * RPMean(pp[ip - sy_p], pp[ip],
                                              rpp[ip - sy_p] * dp[ip - sy_p], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope *
                                      PMean(pp[ip], pp[ip - sy_p], permyp[ip],
                                            permyp[ip - sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip - sy_p], rpp[ip] * dp[ip],
                                               rpp[ip - sy_p] * dp[ip - sy_p])
                                      / viscosity;


                             diff = value - pp[ip];
                             u_new = RPMean(value, pp[ip],
                                            rpp[ip - sy_p] * dp[ip - sy_p], rpp[ip] * dp[ip]);
                             u_new = u_new * z_mult_dat[ip] * ffy * del_x_slope * (permyp[ip] / viscosity)
                                     * 2.0 * (diff / dy);
                           }),
                           FACE(Up,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + sy_p];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip], pp[ip + sy_p],
                                             permyp[ip], permyp[ip + sy_p])
                                     * (diff / dy * del_y_slope)
                                     * RPMean(pp[ip], pp[ip + sy_p],
                                              rpp[ip] * dp[ip], rpp[ip + sy_p] * dp[ip + sy_p])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope
                                      * PMean(pp[ip], pp[ip + sy_p], permyp[ip],
                                              permyp[ip + sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip + sy_p], rpp[ip] * dp[ip],
                                               rpp[ip + sy_p] * dp[ip + sy_p])
                                      / viscosity;


                             diff = pp[ip] - value;
                             u_new = RPMean(pp[ip], value,
                                            rpp[ip] * dp[ip], rpp[ip + sy_p] * dp[ip + sy_p]);
                             u_new = u_new * z_mult_dat[ip] * ffy * del_x_slope * (permyp[ip] / viscosity)
                                     * 2.0 * (diff / dy);
                           }),
                           FACE(Back,
                           {
                             dir = -1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_p]); //RMM

                             lower_cond = pp[ip - sz_p] / sep
                                          - (z_mult_dat[ip - sz_p] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip - sz_p] * gravity *
                                          z_dir_g;

                             upper_cond = pp[ip] / sep + (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip] * gravity *
                                          z_dir_g;

                             diff = (lower_cond - upper_cond);

                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip - sz_p], permzp[ip],
                                               z_mult_dat[ip - sz_p], z_mult_dat[ip])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip - sz_p] * dp[ip - sz_p], rpp[ip] * dp[ip])
                                     / viscosity;

                             sep = dz * z_mult_dat[ip] / 2.0;

                             lower_cond = value / sep - 0.25 * dp[ip] * gravity;
                             upper_cond = pp[ip] / sep + 0.25 * dp[ip] * gravity;
                             diff = (lower_cond - upper_cond);
                             u_new = RPMean(lower_cond, upper_cond,
                                            rpp[ip - sz_p] * dp[ip - sz_p], rpp[ip] * dp[ip]);
                             u_new = u_new * ffz * del_x_slope * del_y_slope *
                                     (permzp[ip] / viscosity)
                                     * 2.0 * diff;
                           }),
                           FACE(Front,
                           {
                             dir = 1;

                             /* Calculate upper face velocity.
                              * @RMM added cos to g term to test terrain-following grid
                              */

                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_p]); //RMM

                             lower_cond = pp[ip] / sep - (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip] * gravity *
                                          z_dir_g;

                             upper_cond = pp[ip + sz_p] / sep
                                          + (z_mult_dat[ip + sz_p] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p]))
                                          * dp[ip + sz_p] * gravity * z_dir_g;

                             diff = (lower_cond - upper_cond);


                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip], permzp[ip + sz_p],
                                               z_mult_dat[ip], z_mult_dat[ip + sz_p])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip] * dp[ip], rpp[ip + sz_p] * dp[ip + sz_p])
                                     / viscosity;

                             sep = dz * z_mult_dat[ip] / 2.0;

                             lower_cond = (pp[ip] / sep) - 0.25 * dp[ip] * gravity * z_dir_g;
                             upper_cond = (value / sep) + 0.25 * dp[ip] * gravity * z_dir_g;

                             diff = lower_cond - upper_cond;
                             u_new = RPMean(lower_cond, upper_cond,
                                            rpp[ip] * dp[ip], rpp[ip + sz_p] * dp[ip + sz_p]);
                             u_new = u_new * ffz * del_x_slope * del_y_slope *
                                     (permzp[ip] / viscosity)
                                     * 2.0 * diff;
                           }),
                           CellFinalize(
                           {
/* Remove the boundary term computed above */
                             fp[ip] -= dt * dir * u_old;

                             /* Add the correct boundary term */
                             fp[ip] += dt * dir * u_new;
                           }),
                           AfterAllCells(DoNothing)
        ); /* End DirichletBC */

      ForPatchCellsPerFace(FluxBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup(
                           {
                             ip = SubvectorEltIndex(p_sub, i, j, k);
                             io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

                             x_dir_g = 0.0;
                             y_dir_g = 0.0;
                             z_dir_g = 1.0;

                             del_x_slope = 1.0;
                             del_y_slope = 1.0;
                           }),
                           FACE(Left,
                           {
                             dir = -1;
                             diff = pp[ip - 1] - pp[ip];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip - 1], pp[ip],
                                             permxp[ip - 1], permxp[ip])
                                     * (diff / dx * del_x_slope)
                                     * RPMean(pp[ip - 1], pp[ip],
                                              rpp[ip - 1] * dp[ip - 1], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip - 1], pp[ip],
                                              permxp[ip - 1], permxp[ip])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip - 1], pp[ip], rpp[ip - 1] * dp[ip - 1],
                                               rpp[ip] * dp[ip])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx;
                           }),
                           FACE(Right,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + 1];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip], pp[ip + 1],
                                             permxp[ip], permxp[ip + 1])
                                     * (diff / dx * del_x_slope)
                                     * RPMean(pp[ip], pp[ip + 1],
                                              rpp[ip] * dp[ip], rpp[ip + 1] * dp[ip + 1])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip], pp[ip + 1],
                                              permxp[ip], permxp[ip + 1])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip], pp[ip + 1], rpp[ip] * dp[ip],
                                               rpp[ip + 1] * dp[ip + 1])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx;
                           }),
                           FACE(Down,
                           {
                             dir = -1;
                             diff = pp[ip - sy_p] - pp[ip];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip - sy_p], pp[ip],
                                             permyp[ip - sy_p], permyp[ip])
                                     * (diff / dy)
                                     * RPMean(pp[ip - sy_p], pp[ip],
                                              rpp[ip - sy_p] * dp[ip - sy_p], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope *
                                      PMean(pp[ip], pp[ip - sy_p], permyp[ip],
                                            permyp[ip - sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip - sy_p], rpp[ip] * dp[ip],
                                               rpp[ip - sy_p] * dp[ip - sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Up,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + sy_p];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip], pp[ip + sy_p],
                                             permyp[ip], permyp[ip + sy_p])
                                     * (diff / dy)
                                     * RPMean(pp[ip], pp[ip + sy_p],
                                              rpp[ip] * dp[ip], rpp[ip + sy_p] * dp[ip + sy_p])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope
                                      * PMean(pp[ip], pp[ip + sy_p], permyp[ip],
                                              permyp[ip + sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip + sy_p], rpp[ip] * dp[ip],
                                               rpp[ip + sy_p] * dp[ip + sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Back,
                           {
                             dir = -1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_p]); //RMM

                             lower_cond = (pp[ip - sz_p] / sep)
                                          - (z_mult_dat[ip - sz_p] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip - sz_p] * gravity *
                                          z_dir_g;

                             upper_cond = (pp[ip] / sep) + (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip] * gravity *
                                          z_dir_g;

                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip - sz_p], permzp[ip],
                                               z_mult_dat[ip - sz_p], z_mult_dat[ip])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip - sz_p] * dp[ip - sz_p], rpp[ip] * dp[ip])
                                     / viscosity;
                             u_new = ffz * del_x_slope * del_y_slope;
                           }),
                           FACE(Front,
                           {
                             dir = 1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_p]); //RMM

                             lower_cond = (pp[ip] / sep) - (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip] * gravity *
                                          z_dir_g;

                             upper_cond = (pp[ip + sz_p] / sep)
                                          + (z_mult_dat[ip + sz_p] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip + sz_p] * gravity *
                                          z_dir_g;

                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip], permzp[ip + sz_p],
                                               z_mult_dat[ip], z_mult_dat[ip + sz_p])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip] * dp[ip], rpp[ip + sz_p] * dp[ip + sz_p])
                                     / viscosity;
                             u_new = ffz * del_x_slope * del_y_slope;
                           }),
                           CellFinalize(
                           {
/* Remove the boundary term computed above */
                             fp[ip] -= dt * dir * u_old;
                             /* Add the correct boundary term */
                             u_new = u_new * bc_patch_values[ival];
                             fp[ip] += dt * dir * u_new;
                           }),
                           AfterAllCells(DoNothing)
        ); /* End FluxBC */

      ForPatchCellsPerFace(OverlandBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(
                           {
                             if (diffusive == 0)
                             {
                               /* Call overlandflow_eval to compute fluxes across the east, west, north, and south faces */
                               PFModuleInvokeType(OverlandFlowEvalInvoke, overlandflow_module,
                                                  (grid, is, bc_struct, ipatch,
                                                   problem_data, pressure, old_pressure,
                                                   ke_, kw_, kn_, ks_, qx_, qy_, CALCFCN));
                             }
                             else
                             {
                               /*  @RMM this is modified to be kinematic wave routing, with a new module for diffusive wave
                                * routing added */
                               double *dummy1;
                               double *dummy2;
                               double *dummy3;
                               double *dummy4;
                               PFModuleInvokeType(OverlandFlowEvalDiffInvoke, overlandflow_module_diff,
                                                  (grid, is, bc_struct, ipatch,
                                                   problem_data, pressure, old_pressure,
                                                   ke_, kw_, kn_, ks_,
                                                   dummy1, dummy2, dummy3, dummy4,
                                                   qx_, qy_, CALCFCN));
                             }
                           }),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup(
                           {
                             ip = SubvectorEltIndex(p_sub, i, j, k);
                             io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

                             x_dir_g = 0.0;
                             y_dir_g = 0.0;
                             z_dir_g = 1.0;

                             del_x_slope = 1.0;
                             del_y_slope = 1.0;
                           }),
                           FACE(Left,
                           {
                             dir = -1;
                             diff = pp[ip - 1] - pp[ip];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip - 1], pp[ip],
                                             permxp[ip - 1], permxp[ip])
                                     * (diff / dx)
                                     * RPMean(pp[ip - 1], pp[ip],
                                              rpp[ip - 1] * dp[ip - 1], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip - 1], pp[ip],
                                              permxp[ip - 1], permxp[ip])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip - 1], pp[ip], rpp[ip - 1] * dp[ip - 1],
                                               rpp[ip] * dp[ip])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx * del_y_slope;
                           }),
                           FACE(Right,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + 1];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip], pp[ip + 1],
                                             permxp[ip], permxp[ip + 1])
                                     * (diff / dx)
                                     * RPMean(pp[ip], pp[ip + 1],
                                              rpp[ip] * dp[ip], rpp[ip + 1] * dp[ip + 1])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip], pp[ip + 1],
                                              permxp[ip], permxp[ip + 1])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip], pp[ip + 1], rpp[ip] * dp[ip],
                                               rpp[ip + 1] * dp[ip + 1])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx * del_y_slope;
                           }),
                           FACE(Down,
                           {
                             dir = -1;
                             diff = pp[ip - sy_p] - pp[ip];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip - sy_p], pp[ip],
                                             permyp[ip - sy_p], permyp[ip])
                                     * (diff / dy)
                                     * RPMean(pp[ip - sy_p], pp[ip],
                                              rpp[ip - sy_p] * dp[ip - sy_p], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope *
                                      PMean(pp[ip], pp[ip - sy_p], permyp[ip],
                                            permyp[ip - sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip - sy_p], rpp[ip] * dp[ip],
                                               rpp[ip - sy_p] * dp[ip - sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Up,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + sy_p];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip], pp[ip + sy_p],
                                             permyp[ip], permyp[ip + sy_p])
                                     * (diff / dy)
                                     * RPMean(pp[ip], pp[ip + sy_p],
                                              rpp[ip] * dp[ip], rpp[ip + sy_p] * dp[ip + sy_p])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope
                                      * PMean(pp[ip], pp[ip + sy_p], permyp[ip],
                                              permyp[ip + sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip + sy_p], rpp[ip] * dp[ip],
                                               rpp[ip + sy_p] * dp[ip + sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Back,
                           {
                             dir = -1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_p]); //RMM
                             //  sep = dz*z_mult_dat[ip];  //RMM

                             lower_cond = (pp[ip - sz_p] / sep)
                                          - (z_mult_dat[ip - sz_p] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip - sz_p] * gravity *
                                          z_dir_g;
                             upper_cond = (pp[ip] / sep) + (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip] * gravity *
                                          z_dir_g;

                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip - sz_p], permzp[ip],
                                               z_mult_dat[ip - sz_p], z_mult_dat[ip])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip - sz_p] * dp[ip - sz_p], rpp[ip] * dp[ip])
                                     / viscosity;
                             u_new = ffz * del_x_slope * del_y_slope;
                           }),
                           FACE(Front,
                           {
                             dir = 1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_p]); //RMM

                             lower_cond = (pp[ip] / sep) - (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip] * gravity *
                                          z_dir_g;
                             upper_cond = (pp[ip + sz_p] / sep)
                                          + (z_mult_dat[ip + sz_p] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip + sz_p] * gravity *
                                          z_dir_g;
                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip], permzp[ip + sz_p],
                                               z_mult_dat[ip], z_mult_dat[ip + sz_p])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip] * dp[ip], rpp[ip + sz_p] * dp[ip + sz_p])
                                     / viscosity;
                             u_new = ffz * del_x_slope * del_y_slope;

                             /* Add overland contribs */
                             io = SubvectorEltIndex(x_sl_sub, i, j, 0);

                             q_overlnd = 0.0;
                             q_overlnd = vol
                                         * (pfmax(pp[ip], 0.0) - pfmax(opp[ip], 0.0)) / dz +
                                         dt * vol * ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy)
                                         / dz + vol * dt / dz * (exp(pfmin(pp[ip], 0.0) * public_xtra->SpinupDampP1) * public_xtra->SpinupDampP2);
                             //NBE

                             if (overlandspinup == 1)
                             {
                               /* add flux loss equal to excess head  that overwrites the prior overland flux */
                               q_overlnd = (vol / dz) * dt
                                           * ((pfmax(pp[ip], 0.0) - 0.0)
                                              + exp(pfmin(pp[ip], 0.0) * public_xtra->SpinupDampP1)
                                              * public_xtra->SpinupDampP2); //@RMM
                             }
                             fp[ip] += q_overlnd;
                           }),
                           CellFinalize(
                           {
/* Remove the boundary term computed above */
                             fp[ip] -= dt * dir * u_old;
                             //add source boundary terms
                             u_new = u_new * bc_patch_values[ival];       //sk: here we go in and implement surface routing!

                             fp[ip] += dt * dir * u_new;
                           }),
                           AfterAllCells(DoNothing)
        ); /* End OverlandBC case */

      ForPatchCellsPerFace(SeepageFaceBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup(
                           {
                             ip = SubvectorEltIndex(p_sub, i, j, k);
                             io = SubmatrixEltIndex(x_ssl_sub, i, j, grid2d_iz);

                             x_dir_g = 0.0;
                             y_dir_g = 0.0;
                             z_dir_g = 1.0;

                             del_x_slope = 1.0;
                             del_y_slope = 1.0;
                             u_old = 0.0;
                           }),
                           FACE(Left, {
                               dir = -1;
                               u_new = z_mult_dat[ip] * ffx * del_y_slope;
                             }),
                           FACE(Right, {
                               dir = 1;
                               u_new = z_mult_dat[ip] * ffx * del_y_slope;
                             }),
                           FACE(Down, {
                               dir = -1;
                               u_new = z_mult_dat[ip] * ffy * del_x_slope;
                             }),
                           FACE(Up, {
                               dir = 1;
                               u_new = z_mult_dat[ip] * ffy * del_x_slope;
                             }),
                           FACE(Back, {
                               dir = -1;
                               u_new = ffz * del_x_slope * del_y_slope;
                             }),
                           FACE(Front, {
                               dir = 1;
                               u_new = ffz * del_x_slope * del_y_slope;

                               io = SubmatrixEltIndex(x_sl_sub, i, j, 0);
/* add flux loss equal to excess head that overwrites the prior overland flux */
                               q_overlnd = (vol / dz) * dt * (pfmax(pp[ip], 0.0) - 0.0); //@RMM

                               fp[ip] += q_overlnd;
                             }),
                           CellFinalize(
                           {
/* Remove the boundary condition computed above */
                             fp[ip] -= dt * dir * u_old;
                             // add source boundary terms
                             u_new = u_new * bc_patch_values[ival];
                             fp[ip] += dt * dir * u_new;
                           }),
                           AfterAllCells(DoNothing)
        ); /* End SeepageFaceBC case */

      ForPatchCellsPerFace(OverlandKinematicBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(
                           {
/*  @RMM this is modified to be kinematic wave routing, with a new module for diffusive wave
 * routing added */
                             double *dummy1;
                             double *dummy2;
                             double *dummy3;
                             double *dummy4;
                             PFModuleInvokeType(OverlandFlowEvalKinInvoke, overlandflow_module_kin,
                                                (grid, is, bc_struct, ipatch, problem_data, pressure,
                                                 ke_, kw_, kn_, ks_,
                                                 dummy1, dummy2, dummy3, dummy4,
                                                 qx_, qy_, CALCFCN));
                           }),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup(
                           {
                             ip = SubvectorEltIndex(p_sub, i, j, k);
                             io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

                             x_dir_g = 0.0;
                             y_dir_g = 0.0;
                             z_dir_g = 1.0;

                             del_x_slope = 1.0;
                             del_y_slope = 1.0;
                           }),
                           FACE(Left,
                           {
                             dir = -1;
                             diff = pp[ip - 1] - pp[ip];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip - 1], pp[ip],
                                             permxp[ip - 1], permxp[ip])
                                     * (diff / dx)
                                     * RPMean(pp[ip - 1], pp[ip],
                                              rpp[ip - 1] * dp[ip - 1], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip - 1], pp[ip],
                                              permxp[ip - 1], permxp[ip])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip - 1], pp[ip], rpp[ip - 1] * dp[ip - 1],
                                               rpp[ip] * dp[ip])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx * del_y_slope;
                           }),
                           FACE(Right,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + 1];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip], pp[ip + 1],
                                             permxp[ip], permxp[ip + 1])
                                     * (diff / dx)
                                     * RPMean(pp[ip], pp[ip + 1],
                                              rpp[ip] * dp[ip], rpp[ip + 1] * dp[ip + 1])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip], pp[ip + 1],
                                              permxp[ip], permxp[ip + 1])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip], pp[ip + 1], rpp[ip] * dp[ip],
                                               rpp[ip + 1] * dp[ip + 1])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx * del_y_slope;
                           }),
                           FACE(Down,
                           {
                             dir = -1;
                             diff = pp[ip - sy_p] - pp[ip];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip - sy_p], pp[ip],
                                             permyp[ip - sy_p], permyp[ip])
                                     * (diff / dy)
                                     * RPMean(pp[ip - sy_p], pp[ip],
                                              rpp[ip - sy_p] * dp[ip - sy_p], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope *
                                      PMean(pp[ip], pp[ip - sy_p], permyp[ip],
                                            permyp[ip - sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip - sy_p], rpp[ip] * dp[ip],
                                               rpp[ip - sy_p] * dp[ip - sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Up,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + sy_p];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip], pp[ip + sy_p],
                                             permyp[ip], permyp[ip + sy_p])
                                     * (diff / dy)
                                     * RPMean(pp[ip], pp[ip + sy_p],
                                              rpp[ip] * dp[ip], rpp[ip + sy_p] * dp[ip + sy_p])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope
                                      * PMean(pp[ip], pp[ip + sy_p], permyp[ip],
                                              permyp[ip + sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip + sy_p], rpp[ip] * dp[ip],
                                               rpp[ip + sy_p] * dp[ip + sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Back,
                           {
                             dir = -1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_p]); //RMM
                             //  sep = dz*z_mult_dat[ip];  //RMM

                             lower_cond = (pp[ip - sz_p] / sep)
                                          - (z_mult_dat[ip - sz_p] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip - sz_p] * gravity *
                                          z_dir_g;
                             upper_cond = (pp[ip] / sep) + (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip] * gravity *
                                          z_dir_g;

                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip - sz_p], permzp[ip],
                                               z_mult_dat[ip - sz_p], z_mult_dat[ip])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip - sz_p] * dp[ip - sz_p], rpp[ip] * dp[ip])
                                     / viscosity;
                             u_new = ffz * del_x_slope * del_y_slope;
                           }),
                           FACE(Front,
                           {
                             dir = 1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_p]); //RMM

                             lower_cond = (pp[ip] / sep) - (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip] * gravity *
                                          z_dir_g;
                             upper_cond = (pp[ip + sz_p] / sep)
                                          + (z_mult_dat[ip + sz_p] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip + sz_p] * gravity *
                                          z_dir_g;
                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip], permzp[ip + sz_p],
                                               z_mult_dat[ip], z_mult_dat[ip + sz_p])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip] * dp[ip], rpp[ip + sz_p] * dp[ip + sz_p])
                                     / viscosity;
                             u_new = ffz * del_x_slope * del_y_slope;

                             io = SubvectorEltIndex(x_sl_sub, i, j, 0);
                             q_overlnd = 0.0;
                             q_overlnd = vol
                                         * (pfmax(pp[ip], 0.0) - pfmax(opp[ip], 0.0)) / dz +
                                         dt * vol * ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy)
                                         / dz;
                             fp[ip] += q_overlnd;
                           }),
                           CellFinalize(
                           {
/* Remove the boundary term computed above */
                             fp[ip] -= dt * dir * u_old;
                             //add source boundary terms
                             u_new = u_new * bc_patch_values[ival];

                             fp[ip] += dt * dir * u_new;
                           }),
                           AfterAllCells(DoNothing)
        );

      ForPatchCellsPerFace(OverlandDiffusiveBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(
                           {
                             /*  @RMM this is a new module for diffusive wave
                              */
                             double *dummy1;
                             double *dummy2;
                             double *dummy3;
                             double *dummy4;
                             PFModuleInvokeType(OverlandFlowEvalDiffInvoke, overlandflow_module_diff,
                                                (grid, is, bc_struct, ipatch,
                                                 problem_data, pressure, old_pressure,
                                                 ke_, kw_, kn_, ks_,
                                                 dummy1, dummy2, dummy3, dummy4,
                                                 qx_, qy_, CALCFCN));
                           }),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup(
                           {
                             ip = SubvectorEltIndex(p_sub, i, j, k);
                             io = SubvectorEltIndex(x_ssl_sub, i, j, grid2d_iz);

                             x_dir_g = 0.0;
                             y_dir_g = 0.0;
                             z_dir_g = 1.0;

                             del_x_slope = 1.0;
                             del_y_slope = 1.0;
                           }),
                           FACE(Left,
                           {
                             dir = -1;
                             diff = pp[ip - 1] - pp[ip];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip - 1], pp[ip],
                                             permxp[ip - 1], permxp[ip])
                                     * (diff / dx)
                                     * RPMean(pp[ip - 1], pp[ip],
                                              rpp[ip - 1] * dp[ip - 1], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip - 1], pp[ip],
                                              permxp[ip - 1], permxp[ip])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip - 1], pp[ip], rpp[ip - 1] * dp[ip - 1],
                                               rpp[ip] * dp[ip])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx * del_y_slope;
                           }),
                           FACE(Right,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + 1];
                             u_old = z_mult_dat[ip] * ffx * del_y_slope
                                     * PMean(pp[ip], pp[ip + 1],
                                             permxp[ip], permxp[ip + 1])
                                     * (diff / dx)
                                     * RPMean(pp[ip], pp[ip + 1],
                                              rpp[ip] * dp[ip], rpp[ip + 1] * dp[ip + 1])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffx * del_y_slope
                                      * PMean(pp[ip], pp[ip + 1],
                                              permxp[ip], permxp[ip + 1])
                                      * (-x_dir_g)
                                      * RPMean(pp[ip], pp[ip + 1], rpp[ip] * dp[ip],
                                               rpp[ip + 1] * dp[ip + 1])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffx * del_y_slope;
                           }),
                           FACE(Down,
                           {
                             dir = -1;
                             diff = pp[ip - sy_p] - pp[ip];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip - sy_p], pp[ip],
                                             permyp[ip - sy_p], permyp[ip])
                                     * (diff / dy)
                                     * RPMean(pp[ip - sy_p], pp[ip],
                                              rpp[ip - sy_p] * dp[ip - sy_p], rpp[ip] * dp[ip])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope *
                                      PMean(pp[ip], pp[ip - sy_p], permyp[ip],
                                            permyp[ip - sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip - sy_p], rpp[ip] * dp[ip],
                                               rpp[ip - sy_p] * dp[ip - sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Up,
                           {
                             dir = 1;
                             diff = pp[ip] - pp[ip + sy_p];
                             u_old = z_mult_dat[ip] * ffy * del_x_slope
                                     * PMean(pp[ip], pp[ip + sy_p],
                                             permyp[ip], permyp[ip + sy_p])
                                     * (diff / dy)
                                     * RPMean(pp[ip], pp[ip + sy_p],
                                              rpp[ip] * dp[ip], rpp[ip + sy_p] * dp[ip + sy_p])
                                     / viscosity;

                             u_old += z_mult_dat[ip] * ffy * del_x_slope
                                      * PMean(pp[ip], pp[ip + sy_p], permyp[ip],
                                              permyp[ip + sy_p])
                                      * (-y_dir_g)
                                      * RPMean(pp[ip], pp[ip + sy_p], rpp[ip] * dp[ip],
                                               rpp[ip + sy_p] * dp[ip + sy_p])
                                      / viscosity;
                             u_new = z_mult_dat[ip] * ffy * del_x_slope;
                           }),
                           FACE(Back,
                           {
                             dir = -1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_p]); //RMM
                             //  sep = dz*z_mult_dat[ip];  //RMM

                             lower_cond = (pp[ip - sz_p] / sep)
                                          - (z_mult_dat[ip - sz_p] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip - sz_p] * gravity *
                                          z_dir_g;
                             upper_cond = (pp[ip] / sep) + (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip - sz_p])) * dp[ip] * gravity *
                                          z_dir_g;

                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip - sz_p], permzp[ip],
                                               z_mult_dat[ip - sz_p], z_mult_dat[ip])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip - sz_p] * dp[ip - sz_p], rpp[ip] * dp[ip])
                                     / viscosity;
                           }),
                           FACE(Front,
                           {
                             dir = 1;
                             sep = dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_p]); //RMM

                             lower_cond = (pp[ip] / sep) - (z_mult_dat[ip] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip] * gravity *
                                          z_dir_g;
                             upper_cond = (pp[ip + sz_p] / sep)
                                          + (z_mult_dat[ip + sz_p] / (z_mult_dat[ip] + z_mult_dat[ip + sz_p])) * dp[ip + sz_p] * gravity *
                                          z_dir_g;
                             diff = lower_cond - upper_cond;
                             u_old = ffz * del_x_slope * del_y_slope
                                     * PMeanDZ(permzp[ip], permzp[ip + sz_p],
                                               z_mult_dat[ip], z_mult_dat[ip + sz_p])
                                     * diff
                                     * RPMean(lower_cond, upper_cond,
                                              rpp[ip] * dp[ip], rpp[ip + sz_p] * dp[ip + sz_p])
                                     / viscosity;

                             io = SubvectorEltIndex(x_sl_sub, i, j, 0);
                             q_overlnd = 0.0;
                             q_overlnd = vol
                                         * (pfmax(pp[ip], 0.0) - pfmax(opp[ip], 0.0)) / dz +
                                         dt * vol * ((ke_[io] - kw_[io]) / dx + (kn_[io] - ks_[io]) / dy)
                                         / dz;

                             fp[ip] += q_overlnd;
                           }),
                           CellFinalize(
                           {
                             /* Remove the boundary term computed above */
                             fp[ip] -= dt * dir * u_old;
                             //add source boundary terms
                             u_new = u_new * bc_patch_values[ival];       //sk: here we go in and implement surface routing!
                             fp[ip] += dt * dir * u_new;
                           }),
                           AfterAllCells(DoNothing)
        ); /* End OverlandDiffusiveBC case */
    }          /* End ipatch loop */
  }            /* End subgrid loop */

  /*
   * Reset values inserted for the DirichletBC back to the decoupled
   * problem used in the inactive cells.
   *
   * See comments above on why this is needed.
   */
  ForSubgridI(is, GridSubgrids(grid))
  {
    subgrid = GridSubgrid(grid, is);

    p_sub = VectorSubvector(pressure, is);
    f_sub = VectorSubvector(fval, is);

    nx_p = SubvectorNX(p_sub);
    ny_p = SubvectorNY(p_sub);
    nz_p = SubvectorNZ(p_sub);

    sy_p = nx_p;
    sz_p = ny_p * nx_p;

    pp = SubvectorData(p_sub);
    fp = SubvectorData(f_sub);

    int pp_idx = 0;
    ForBCStructNumPatches(ipatch, bc_struct)
    {
      bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);

      ForPatchCellsPerFace(DirichletBC,
                           InParallel, NO_LOCALS,
                           BeforeAllCells(DoNothing),
                           LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                           CellSetup({
                               pp_idx = 0;
                               ip = SubvectorEltIndex(p_sub, i, j, k);
                               value = bc_patch_values[ival];
                             }),
                           FACE(Left,  { pp_idx = ip - 1; }),
                           FACE(Right, { pp_idx = ip + 1; }),
                           FACE(Down,  { pp_idx = ip - sy_p; }),
                           FACE(Up,    { pp_idx = ip + sy_p; }),
                           FACE(Back,  { pp_idx = ip - sz_p; }),
                           FACE(Front, { pp_idx = ip + sz_p; }),
                           CellFinalize({
                               pp[pp_idx] = -FLT_MAX;
                               fp[pp_idx] = 0.0;
                             }),
                           AfterAllCells(DoNothing)
        );
    }          /* End ipatch loop */
  }            /* End subgrid loop */


  } // End Parallel region

  FreeBCStruct(bc_struct);


  PFModuleInvokeType(RichardsBCInternalInvoke, bc_internal, (problem, problem_data, fval, NULL,
                                                             time, pressure, CALCFCN));

  EndTiming(public_xtra->time_index);

  FreeVector(KW);
  FreeVector(KE);
  FreeVector(KN);
  FreeVector(KS);
  FreeVector(qx);
  FreeVector(qy);

  FreeVector(u_right_vec);
  FreeVector(u_front_vec);
  FreeVector(u_upper_vec);

  return;
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule    *NlFunctionEvalInitInstanceXtra(Problem *problem,
                                            Grid *   grid,
                                            double * temp_data)

{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra;

  (void)grid;
  (void)temp_data;

  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (problem != NULL)
  {
    (instance_xtra->problem) = problem;
  }

  if (PFModuleInstanceXtra(this_module) == NULL)
  {
    (instance_xtra->density_module) =
      PFModuleNewInstance(ProblemPhaseDensity(problem), ());
    (instance_xtra->saturation_module) =
      PFModuleNewInstanceType(SaturationInitInstanceXtraInvoke,
                              ProblemSaturation(problem), (NULL, NULL));
    (instance_xtra->rel_perm_module) =
      PFModuleNewInstanceType(PhaseRelPermInitInstanceXtraInvoke,
                              ProblemPhaseRelPerm(problem), (NULL, NULL));
    (instance_xtra->phase_source) =
      PFModuleNewInstance(ProblemPhaseSource(problem), ());
    (instance_xtra->bc_pressure) =
      PFModuleNewInstanceType(BCPressurePackageInitInstanceXtraInvoke,
                              ProblemBCPressure(problem), (problem));
    (instance_xtra->bc_internal) =
      PFModuleNewInstance(ProblemBCInternal(problem), ());
    (instance_xtra->overlandflow_module) =
      PFModuleNewInstance(ProblemOverlandFlowEval(problem), ());     //DOK
    (instance_xtra->overlandflow_module_diff) =
      PFModuleNewInstance(ProblemOverlandFlowEvalDiff(problem), ());   //@RMM
    (instance_xtra->overlandflow_module_kin) =
      PFModuleNewInstance(ProblemOverlandFlowEvalKin(problem), ());
  }
  else
  {
    PFModuleReNewInstance((instance_xtra->density_module), ());
    PFModuleReNewInstanceType(SaturationInitInstanceXtraInvoke,
                              (instance_xtra->saturation_module),
                              (NULL, NULL));
    PFModuleReNewInstanceType(PhaseRelPermInitInstanceXtraInvoke,
                              (instance_xtra->rel_perm_module),
                              (NULL, NULL));
    PFModuleReNewInstanceType(BCPressurePackageInitInstanceXtraInvoke,
                              (instance_xtra->phase_source), (NULL));
    PFModuleReNewInstanceType(BCPressurePackageInitInstanceXtraInvoke,
                              (instance_xtra->bc_pressure), (problem));
    PFModuleReNewInstance((instance_xtra->bc_internal), ());
    PFModuleReNewInstance((instance_xtra->overlandflow_module), ());     //DOK
    PFModuleReNewInstance((instance_xtra->overlandflow_module_diff), ());      //@RMM
    PFModuleReNewInstance((instance_xtra->overlandflow_module_kin), ());
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  NlFunctionEvalFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    PFModuleFreeInstance(instance_xtra->density_module);
    PFModuleFreeInstance(instance_xtra->saturation_module);
    PFModuleFreeInstance(instance_xtra->rel_perm_module);
    PFModuleFreeInstance(instance_xtra->phase_source);
    PFModuleFreeInstance(instance_xtra->bc_pressure);
    PFModuleFreeInstance(instance_xtra->bc_internal);
    PFModuleFreeInstance(instance_xtra->overlandflow_module);     //DOK
    PFModuleFreeInstance(instance_xtra->overlandflow_module_diff);      //@RMM
    PFModuleFreeInstance(instance_xtra->overlandflow_module_kin);

    tfree(instance_xtra);
  }
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule   *NlFunctionEvalNewPublicXtra(char *name)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;
  char key[IDB_MAX_KEY_LEN];
  char *switch_name;
  int switch_value;
  NameArray switch_na;
  NameArray upwind_switch_na;


  public_xtra = ctalloc(PublicXtra, 1);

/* These parameters dampen the transition/switching into overland flow to speedup
 * the spinup process. */
  sprintf(key, "OverlandSpinupDampP1");
  public_xtra->SpinupDampP1 = GetDoubleDefault(key, 0.0);
  sprintf(key, "OverlandSpinupDampP2");
  public_xtra->SpinupDampP2 = GetDoubleDefault(key, 0.0);    //NBE

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

  (public_xtra->time_index) = RegisterTiming("NL_F_Eval");

  PFModulePublicXtra(this_module) = public_xtra;

  return this_module;
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalFreePublicXtra
 *--------------------------------------------------------------------------*/

void  NlFunctionEvalFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);


  if (public_xtra)
  {
    tfree(public_xtra);
  }
}


/*--------------------------------------------------------------------------
 * NlFunctionEvalSizeOfTempData
 *--------------------------------------------------------------------------*/

int  NlFunctionEvalSizeOfTempData()
{
  return 0;
}

#ifdef USING_PARALLEL
} // Extern C
#endif
